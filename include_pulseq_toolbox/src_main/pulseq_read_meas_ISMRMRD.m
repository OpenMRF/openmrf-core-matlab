function [rawdata, study_info, PULSEQ] = pulseq_read_meas_ISMRMRD(study)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 25.08.2026

% input: only study, which is the path to an ISMRMRD .h5 file
%
% Requires the ISMRMRD MATLAB package (+ismrmrd).
%
% ISMRMRD contents used:
%   data:          MRI rawdata and acquisition headers
%   protocolName:  Pulseq sequence name (YYMMDD_HHMM_user_sequence)
%   H1resonanceFrequency_Hz: larmor frequency [Hz]
%   acquisition_time_stamp: ADC time stamps (Siemens clock ticks)

%% check for file existence and emit clear exception upon missing file
study = strrep(study, '\', '/');
[temp_path, temp_name, temp_ext] = fileparts(study);
if isempty(temp_ext)
    study = [temp_path '/' temp_name '.h5'];
end
clear temp_path temp_name temp_ext;
checkexists(study, file_qualifier='scanner raw data');

%% read ISMRMRD data
try
    dset = ismrmrd.Dataset(study);
catch
    error('ISMRMRD data could not be opened. Make sure that the ISMRMRD MATLAB package is on the MATLAB path!')
end

acq = dset.readAcquisition();

NR = numel(acq.data);
if NR==0
    error('no ISMRMRD acquisitions found!')
end

NRead  = double(acq.head.number_of_samples(1));
NCoils = double(acq.head.active_channels(1));

if any(double(acq.head.number_of_samples)~=NRead)
    error('number of ADC samples is not constant!')
end
if any(double(acq.head.active_channels)~=NCoils)
    error('number of active receiver channels is not constant!')
end

rawdata = zeros(NCoils, NR, NRead, 'single');

for i = 1:NR
    temp = acq.data{i};

    if isequal(size(temp), [NRead NCoils])
        temp = temp.';
    elseif ~isequal(size(temp), [NCoils NRead])
        error('unexpected rawdata dimensions in acquisition %d!', i)
    end

    rawdata(:,i,:) = reshape(single(temp), NCoils, 1, NRead);
end
clear temp;

%% read ISMRMRD XML header
xml_string = dset.readxml();
xml_string = char(xml_string(:).');
xml_string(xml_string==0) = [];

external_name = xml_get(xml_string, 'protocolName');
if isempty(external_name)
    error('protocolName not found in ISMRMRD XML header!')
end

f0 = xml_get(xml_string, 'H1resonanceFrequency_Hz');
if isempty(f0)
    f0 = [];
else
    f0 = str2double(f0);
end

meas_name  = xml_get(xml_string, 'protocolName');
meas_date  = xml_get(xml_string, 'studyDate');
meas_clock = xml_get(xml_string, 'studyTime');

% Siemens acquisition_time_stamp: 2.5 ms clock ticks
time_stamps = double(acq.head.acquisition_time_stamp(:)) * 2.5e-3;

%% import backup of pulseq workspace
[~, external_name] = fileparts(external_name);
[~,temp]      = sort(external_name=='_', 'descend');
seq_id        = external_name(1:11);
scan_id       = seq_id;
scan_id(7)    = [];
scan_id       = int64(str2num(scan_id));
pulseq_user   = external_name( 13 : temp(3)-1 );
seq_name      = external_name( temp(3)+1 : end );
if scan_id==0
    error('scan ID not found!')
end
clear temp;

[~, ~, pulseq_path] = pulseq_get_user_definitions([], 0);
backup_path         = [pulseq_path '/Pulseq_Workspace/' pulseq_user '/' seq_id(1:6) '/' seq_id];
mat_file            = dir(fullfile(backup_path, '*.mat'));
if isempty(mat_file)
    error('no pulseq backup found!')
end
backup_path         = [backup_path '/' mat_file(1).name];
checkexists(backup_path, file_qualifier='pulseq backup');

try
    load(backup_path);
catch
    error('no pulseq backup found!')
end

% old version: rename pulseq_workspace -> PULSEQ
if exist('pulseq_workspace', 'var')
    PULSEQ = pulseq_workspace;
    clear pulseq_workspace;
end

% ISMRMRD conversion does not provide an independent Pulseq MD5 hash
md5_hash = [];

%% output
[study_path, study_name, ext] = fileparts(study);
study_info.hdr                = xml_string;
study_info.seq_name           = seq_name;
study_info.seq_id             = seq_id;
study_info.scan_id            = scan_id;
study_info.md5_hash           = md5_hash;
study_info.pulseq_user        = pulseq_user;
study_info.study_name         = [study_name ext];
study_info.study_path         = study_path;
study_info.backup_file        = backup_path;
study_info.external_name      = external_name;
study_info.meas_name          = meas_name;
study_info.meas_date          = meas_date;
study_info.meas_clock         = meas_clock;
study_info.f0                 = f0;
study_info.time_stamps        = time_stamps(:);
study_info.soft_delays        = [];

end


function value = xml_get(xml_string, tag)

token = regexp(xml_string, ['<' tag '>(.*?)</' tag '>'], 'tokens', 'once');

if isempty(token)
    value = [];
else
    value = token{1};
end

end
