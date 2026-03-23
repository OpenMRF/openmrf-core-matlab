function [rawdata, study_info, PULSEQ] = pulseq_read_meas_UI(study)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 23.03.2026

% input: only study, which is the path to a .mat file containing United Imaging rawdata

% .mat file contents:
%   rawdata:       NCoils x NR x NRead (complex single or double)
%   external_name: the name of the .seq file (e.g. '260322_0426_mgram_mrf.seq' -> YYMMDD_HHMM_user_sequence)
%   md5_hash:      the hash at the end of the .seq file (e.g. '231c7e03ccb5985953a92e213d2a3f77')
%   time_stamps:   NR x 1, adc time stamps [s]
%   f0:            larmor frequency [Hz]
%   soft_delay:    N x 1, [s] list with the soft delay set in the UI interface

%% check for file existence and emit clear exception upon missing file
study = strrep(study, '\', '/');
[temp_path, temp_name, temp_ext] = fileparts(study);
if isempty(temp_ext)
    study = [temp_path '/' temp_name '.mat'];
end
clear temp_path temp_name temp_ext;
checkexists(study, file_qualifier='scanner raw data');

%% check existence of varibales
load(study);

if ~exist('rawdata', 'var')
    error(['rawdata does not exist!']);
end
if ~exist('external_name', 'var')
    error(['external_name does not exist!']);
end
if ~exist('md5_hash', 'var')
    error(['md5_hash does not exist!']);
end
if ~exist('time_stamps', 'var')
    time_stamps = [];
end
if ~exist('f0', 'var')
    f0 = [];
end
if ~exist('soft_delay', 'var')
    soft_delay = [];
end

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
    error('scan ID not found!');
end
clear temp;
[~, ~, pulseq_path] = pulseq_get_user_definitions([], 0);
backup_path         = [pulseq_path '/Pulseq_Workspace/' pulseq_user '/' seq_id(1:6) '/' seq_id];
mat_file            = dir(fullfile(backup_path, '*.mat'));
backup_path         = [backup_path '/' mat_file.name];
checkexists(backup_path, file_qualifier='pulseq backup');
try
    load(backup_path);
catch
    error('no pulseq backup found!')
end

% consistency check: compare md5 hash
if isfield(PULSEQ, 'md5_hash')
    if ~strcmp(md5_hash, PULSEQ.md5_hash)
        warning(' >>>>> MD5 Hashs do not agree! <<<<<')
    end
end

%% output
[study_path, study_name, ext] = fileparts(study);
study_info.hdr                = []; % ask United Imaging
study_info.seq_name           = seq_name;
study_info.seq_id             = seq_id;
study_info.scan_id            = scan_id;
study_info.md5_hash           = md5_hash;
study_info.pulseq_user        = pulseq_user;
study_info.study_name         = [study_name ext];
study_info.study_path         = study_path;
study_info.backup_file        = backup_path;
study_info.external_name      = external_name;
study_info.meas_name          = []; % ask United Imaging
study_info.meas_date          = []; % ask United Imaging
study_info.meas_clock         = []; % ask United Imaging
study_info.f0                 = f0;
study_info.time_stamps        = time_stamps;

end
