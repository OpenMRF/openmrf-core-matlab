function [SEQ, SIM] = MRF_read_seq_file(seq_file, f0, time_stamps, soft_delays, flag_kz, echo_mode, adc_npad, dwell_time, dt, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
% comment regarding pulseq version: with v1.5.1 it is possible to define rotation objects
% the current function ignores rotation objects; this has no effect for
% most sequences but might be an issue to simulate some complicated 3d
% readouts. maybe I will update this in the future...

% ----- input: -----
% seq_file:        path of the .seq file
% f0:              larmor frequency [Hz] for off-resonant pulses or adcs
% time_stamps:     nadc x 1 array of adc time stamps [s]
% soft_delays:     n x 1 array of input delays (e.g. for RR tuning) [s]
% flag_kz:         0 -> no kz paritions; 1 -> search for START/STOP labels and remove unneccesary kz partitions
% echo_mode:       'spiral_out' 'center' 'auto'
% adc_npad:        can be 2x1 or 1x1 and defines the start and end of adc padding
% dwell_time:      dwell time of adc [s]
% dt:              target raster time for the simulation (default 1us)
% flag_plot:       0 -> off; 1 -> calc full waveform arrays; 2 -> plot

% ----- output 1: -----
% SEQ: struct, which contains the full .seq information
% .FILE:        text of the .seq file
% .DEFINITIONS: definitions of the .seq header
% .BLOCKS:      blocks of the .seq file
% .RF:          RF section of the .seq file
% .GRADIENTS:   gradient section of the .seq file
% .TRAP:        trapezodial gradient section of the .seq file
% .ADC:         ADC section of the .seq file
% .EXTENSIONS:  extension list of the .seq file
% .EXT_SPECS:   extension specifications of the .seq file
% .SHAPES:      RF or gradient shapes of the .seq file
% .WAVES:       structured object which contains all waveforms for all blocks
% .FULL:        1d arrays for visualization of sequnce diagrams

% ----- output 2: -----
% SIM: struct, which contains compressed 1d arrays for bloch and epg simulations
%  .ID:  compressed ID array for selecting the correct Bloch operator
%  .RF:  compressed complex RF waveforms [rad/s]
%  .GX:  compressed x gradient moments [rad/m] or waveforms [rad/s/m]
%  .GY:  compressed y gradient moments [rad/m] or waveforms [rad/s/m]
%  .GZ:  compressed z gradient moments [rad/m] or waveforms [rad/s/m]
%  .DB:  compressed b-value of the diffusion encoding tensor [s/mm^2]
%  .DT:  compressed variable time steps [s]
%  .PHI: adc reciever phases [rad]

%% adjust format of seq_file
strrep(seq_file, '\', '/');
if ispc() && strcmp(seq_file(1), '/')
    seq_file = [seq_file(2) ':' seq_file(3:end)];
end
if isunix() && strcmp(seq_file(2), ':')
    seq_file = ['/' seq_file(1) seq_file(3:end)];
end
[~, temp_name] = fileparts(seq_file);
fprintf('\n');
fprintf(['reading pulseq sequence: ' temp_name '.seq \n']);

%% set defaults
if nargin<2
    f0 = [];
end
if nargin<3
    time_stamps = [];
end
if nargin<4
    soft_delays = [];
end
if nargin<5
    flag_kz = [];
end
if nargin<6
    echo_mode = [];
end
if nargin<7
    adc_npad = [];
end
if nargin<8
    dwell_time = [];
end
if nargin<9
    dt = [];
end
if nargin<10
    flag_plot = [];
end
if isempty(f0)
    f0 = 0;
end
if isempty(echo_mode)
    echo_mode = 'spiral_out';
end
if isempty(dt)
    dt = 1e-6;
end
if isempty(flag_plot)
    flag_plot = 0;
end

%% store initial flags and params in SEQ struct
SEQ    = struct();
SEQ.f0 = f0;
SEQ.dt = dt;
SEQ.echo_mode = echo_mode;

%% read main sections of .seq file
tic;
SEQ.FILE        = readlines(seq_file);
SEQ.DEFINITIONS = read_seq_file_definitions(SEQ.FILE);
SEQ.BLOCKS      = read_seq_file_section(SEQ.FILE, 8,  '[BLOCKS]');
SEQ.RF          = read_seq_file_section(SEQ.FILE, 12, '[RF]');
SEQ.GRADIENTS   = read_seq_file_section(SEQ.FILE, 7,  '[GRADIENTS]');
SEQ.TRAP        = read_seq_file_section(SEQ.FILE, 6,  '[TRAP]');
SEQ.ADC         = read_seq_file_section(SEQ.FILE, 9,  '[ADC]');
SEQ.EXTENSIONS  = read_seq_file_section(SEQ.FILE, 4,  '[EXTENSIONS]');
SEQ.EXT_SPECS   = read_seq_file_ext_specs(SEQ.FILE, SEQ.BLOCKS);
SEQ.SHAPES      = read_seq_file_shapes(SEQ.FILE);

%% read and apply extension specifications
SEQ = corr_soft_delays(SEQ, soft_delays);        % correct timings in case of soft delays
SEQ = corr_trigger_delays(SEQ, time_stamps);     % correct timings in case of cardiac trigger
SEQ = filter_kz_partitions(SEQ, flag_kz);        % delete unnecessary blocks in case of kz partitions
SEQ = correct_GE_timings(SEQ);                   % add a 117us delay for each GE TRID label
SEQ.BLOCKS(find(SEQ.BLOCKS(:,2)==0),:) = [];     % delete blocks which only contain labels
SEQ = delete_unnecessary_blocks(SEQ);            % delete blocks before first RF pulse and after last adc

%% calculate waveforms for individual sequence blocks
SEQ = calc_block_waveforms(SEQ);

%% find echo positions within ADCs -> used for signal simulation
if ~isempty(adc_npad) && isempty(dwell_time)
    error('enter dwell_time to correct adc padding');
end
adc_npad = round(adc_npad * dwell_time / dt);
SEQ = find_echo_positions(SEQ, echo_mode, adc_npad);

%% calculate full waveform arrays
if flag_plot>0
    SEQ = calc_full_waveforms(SEQ, dt, flag_plot);
end

%% reduce memory of SEQ_WAVES struct
SEQ = reduce_seq_wave_size(SEQ);

%% init operator list for brute force bloch simulation
temp_t = toc;
fprintf(['   ' num2str(temp_t, '%.1f') 's ... complete! \n']);
tic;
fprintf('initialize operator list \n');
if nargin>1
    SIM = init_brute_force_bloch_simulation(SEQ);
    [~, SIM.seq_name] = fileparts(seq_file);
end
temp_t = toc;
fprintf(['   ' num2str(temp_t, '%.1f') 's ... complete! \n']);
fprintf('\n');

end

% ----------------------------------------------------------------------
% ---------------------------- SUB ROUTINES ----------------------------
% ----------------------------------------------------------------------

%% -------------------- read seq file definitions --------------------
function seq_definitions = read_seq_file_definitions(seq_file)

    seq_definitions.name = ''; 
    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'Name', 4)
            temp = strsplit(seq_file(j));
            seq_definitions.name = temp{2};
            break;
        end
    end
    
    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'Hash', 4)
            temp = strsplit(seq_file(j));
            seq_definitions.hash = temp{2};
            break;
        end
    end

    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'major', 5)
            temp = strsplit(seq_file(j));
            seq_definitions.version_major = str2num(temp{2});
            break;
        end
    end

    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'minor', 5)
            temp = strsplit(seq_file(j));
            seq_definitions.version_minor = str2num(temp{2});
            break;
        end
    end

    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'BlockDurationRaster', 19)
            temp = strsplit(seq_file(j));
            seq_definitions.blockRasterTime = str2num(temp{2});
            break;
        end
    end
    
    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'GradientRasterTime', 18)
            temp = strsplit(seq_file(j));
            seq_definitions.gradRasterTime = str2num(temp{2});
            break;
        end
    end
    
    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'RadiofrequencyRasterTime', 24)
            temp = strsplit(seq_file(j));
            seq_definitions.rfRasterTime = str2num(temp{2});
            break;
        end
    end
    
    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'AdcRasterTime', 13)
            temp = strsplit(seq_file(j));
            seq_definitions.adcRasterTime = str2num(temp{2});
            break;
        end
    end
    
    for j = 1:numel(seq_file)
        if strncmpi(seq_file(j), 'FOV', 3)
            temp = strsplit(seq_file(j));
            seq_definitions.FOV(1) = str2num(temp{2});
            seq_definitions.FOV(2) = str2num(temp{3});
            seq_definitions.FOV(3) = str2num(temp{4});
            break;
        end
    end
    clear j temp;

    if seq_definitions.version_major*10 + seq_definitions.version_minor < 15
        error('compatible only with Pulseq v1.5 or newer!');
    end
    
end

%% -------------------- read seq file sections --------------------
function section = read_seq_file_section(seq_file, n, key)

    % find start of section
    for ind_start = 1:numel(seq_file)
        if strcmp(seq_file(ind_start), key)
            ind_start = ind_start + 1;
            break;
        end
    end

    % find end of section
    for ind_end = ind_start:numel(seq_file)
        if strcmp(seq_file(ind_end), '')
            ind_end = ind_end - 1;
            break;
        end
    end
   
    % read section lines and convert string to numbers
    n_lines = ind_end - ind_start + 1;
    section = zeros(n_lines, n);    
    for j = 1 : n_lines
        if ~strcmp(key, '[RF]')
            section(j,:) = str2num(seq_file(j+ind_start-1));
        else
            temp_line = convertStringsToChars(seq_file(j+ind_start-1));
            switch temp_line(end)
                case 'e'
                    temp_line(end) = '1'; % excitation
                case 'r'
                    temp_line(end) = '2'; % refocusing
                case 'i'
                    temp_line(end) = '3'; % inversion
                case 's'
                    temp_line(end) = '4'; % saturation
                case 'p'
                    temp_line(end) = '5'; % preparation
                otherwise
                    temp_line(end) = '6'; % undefined
            end                
            temp_line = convertCharsToStrings(temp_line);
            section(j,:) = str2num(temp_line);
        end
    end

end

%% -------------------- read seq file extensions specification --------------------
function ext_specs = read_seq_file_ext_specs(seq_file, seq_blocks)

    ext_specs = [];
    if sum(seq_blocks(:,end)) == 0
        return;    
    end

    % find start of section
    for ind_start = 1:numel(seq_file)
        if strcmp(seq_file(ind_start), '[EXTENSIONS]')
            break;
        end
    end
    for ind_start = ind_start:numel(seq_file)
        if strcmp(seq_file(ind_start), '')
            ind_start = ind_start + 1;
            break;
        end
    end
    
    % find end of section
    for ind_end = ind_start:numel(seq_file)
        if strcmp(seq_file(ind_end), '[SHAPES]')
            ind_end = ind_end - 3;
            break;
        end
    end
    
    % read extension specs section
    section = seq_file(ind_start:ind_end);
    clear ind_start ind_end;

    % search for different extension specifications
    ext_specs = struct();
    for j=1:numel(section)
        temp_line = convertStringsToChars(section(j));
        if numel(temp_line)>9
        if strcmp(temp_line(1:9), 'extension')
            ind_start = j;
            ind_end   = numel(section);
            for k = ind_start:numel(section)
                temp_line = convertStringsToChars(section(k));
                if strcmp(temp_line, '')
                    ind_end = k-1;
                    break;
                end
            end
            temp_line = convertStringsToChars(section(ind_start));
            temp_type = str2num(temp_line(end));
            temp_vals = section(ind_start+1:ind_end);
            ext_specs(temp_type,1).type = temp_line(11:end-2);
            ext_specs(temp_type,1).vals = temp_vals;
        end
        end
    end

end

%% -------------------- read seq file shapes --------------------
function seq_shapes = read_seq_file_shapes(seq_file)

    for ind_start = 1:numel(seq_file)
        if strcmp(seq_file(ind_start), '[SHAPES]')
            ind_start = ind_start + 1;
            break;
        end
    end
    for ind_end = 1:numel(seq_file)
        if strcmp(seq_file(ind_end), '[SIGNATURE]')
            ind_end = ind_end + 1;
            break;
        end
    end
    ind_end = ind_end - 2;
    
    shapes_all = seq_file(ind_start:ind_end);
    
    count = 0;
    for j = 1 : numel(shapes_all)
        if strncmpi(shapes_all(j), 'shape_id', 8)
            count = count + 1;
            ind_all(count) = j;
        end
    end
    ind_all(end+1) = numel(shapes_all);
    
    seq_shapes_comp = struct();
    for j = 1 : numel(ind_all)-1
        temp = shapes_all(ind_all(j)+2 : ind_all(j+1)-2);
        for k = 1 : numel(temp)
            seq_shapes_comp(j,1).data(k,1) = str2num(temp(k)); 
        end
        temp = strsplit(shapes_all(ind_all(j)+1));
        temp = temp{1,2};
        seq_shapes_comp(j,1).num_samples = str2num(temp);
    end
    
    seq_shapes = struct();
    for j = 1 : numel(seq_shapes_comp)
        temp = seq_shapes_comp(j,1);
        temp = decompressShape(temp);
        if size(temp,2)>size(temp,1)
            temp = temp';
        end
        seq_shapes(j,1).data = temp;
    end

end

%% -------------------- decompress pulseq shapes --------------------
function w = decompressShape(shape, forceDecompression)
% copied from Pulseq Git: version 20241115
%decompressShape Decompress a gradient or pulse shape.
%   w=decompressShape(shape) Decompress the shape compressed with a run-length
%   compression scheme on the derivative. The given shape is structure with
%   the following fields:
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also compressShape


dataPack = shape.data;
dataPackLen = length(dataPack);
numSamples=shape.num_samples;

if nargin<2
    forceDecompression=false;
end

if ~forceDecompression && numSamples==dataPackLen
    % uncompressed shape
    w=dataPack';
    return;
end

w= zeros(1, numSamples) ;                 % pre-allocate the result matrix
                                          % dimensons: (1,length of the data set)
                                       
% decompression starts here

dataPackDiff = dataPack(2:end) - dataPack(1:end-1);

% when dataPackDiff == 0 the subsequent samples are equal ==> marker for
% repeats (run-length encoding)
dataPackMarkers=find(dataPackDiff==0.0);

countPack= 1;                                               % counter 1: points to the current compressed sample
countUnpack= 1;                                             % counter 2: points to the current uncompressed sample

for i=1:length(dataPackMarkers)
    nextPack=dataPackMarkers(i); % careful, this index may have "false positives" , e.g. if the value 3 repeats 3 times, then we will have 3 3 3
    currUnpackSamples=nextPack-countPack;
    if currUnpackSamples < 0 % this rejects false positives
        continue;        
    elseif currUnpackSamples > 0 % do we have an unpacked block to copy?
        w(countUnpack:(countUnpack+currUnpackSamples-1)) = dataPack(countPack:(nextPack-1));
        countPack = countPack + currUnpackSamples;
        countUnpack = countUnpack + currUnpackSamples;
    end
    % now comes the packed/repeated section
    rep=dataPack(countPack+2)+2;
    w(countUnpack:(countUnpack+rep-1))=dataPack(countPack);
    countPack= countPack + 3;
    countUnpack= countUnpack + rep;
end

% samples left?
if (countPack<=dataPackLen)
    assert(dataPackLen-countPack==numSamples-countUnpack);
    % copy the rest of the shape, it is unpacked
    w(countUnpack:end)= dataPack(countPack:end);
end

w = cumsum(w);
w=w(:);

end

%% -------------------- correct timings in case of soft delays --------------------
function SEQ = corr_soft_delays(SEQ, soft_delays)
    if ~isempty(soft_delays)
        for j = 1:size(SEQ.BLOCKS,1)
            if SEQ.BLOCKS(j,8)>0
                temp_ext_id   = SEQ.BLOCKS(j, 8);
                temp_ext_type = SEQ.EXTENSIONS(temp_ext_id, 2);
                temp_ext_no   = SEQ.EXTENSIONS(temp_ext_id, 3);
                if strcmp(SEQ.EXT_SPECS(temp_ext_type).type, 'DELAYS')
                    temp_ext_vals    = SEQ.EXT_SPECS(temp_ext_type).vals(temp_ext_no);
                    temp_ext_vals    = strsplit(temp_ext_vals, ' ');
                    temp_input       = soft_delays(temp_ext_no); % input [s]
                    temp_offset      = str2num(temp_ext_vals(3))*1e-6; % offset [s]
                    temp_factor      = str2num(temp_ext_vals(4)); % factor
                    temp_dur         = temp_offset + temp_input / temp_factor; % block_duration = offset + input / factor
                    temp_dur         = round(temp_dur/SEQ.DEFINITIONS.blockRasterTime);
                    SEQ.BLOCKS(j, 2) = temp_dur;
                end
            end
        end
    end
end

%% -------------------- correct timings in case of trigger inputs --------------------
function SEQ = corr_trigger_delays(SEQ, adc_time_stamps)
    dt = SEQ.DEFINITIONS.blockRasterTime;    
    if ~isempty(adc_time_stamps)
        % find all blocks with input triggers
        n_trig = 0;
        for j = 1:size(SEQ.BLOCKS,1)
            if SEQ.BLOCKS(j,8)>0
                temp_ext_id   = SEQ.BLOCKS(j, 8);
                temp_ext_type = SEQ.EXTENSIONS(temp_ext_id, 2);
                temp_ext_no   = SEQ.EXTENSIONS(temp_ext_id, 3);
                if strcmp(SEQ.EXT_SPECS(temp_ext_type).type, 'TRIGGERS')
                    temp_ext_vals = SEQ.EXT_SPECS(temp_ext_type).vals(temp_ext_no);
                    temp_ext_vals = str2num(temp_ext_vals);
                    if temp_ext_vals(2) == 2
                        n_trig = n_trig + 1;
                        ind_trig(n_trig,1) = j;
                    end
                end
                clear temp_ext_id temp_ext_type temp_ext_no temp_ext_vals;
            end
        end
        % calculate unknown trigger delays
        if n_trig>0
            for j=1:n_trig
                temp_adcs = SEQ.BLOCKS(:,7)>0;
                temp_adcs_pre = temp_adcs;
                temp_adcs_pre(ind_trig(j):end) = 0;
                temp_adcs_post = temp_adcs;
                temp_adcs_post(1:ind_trig(j)) = 0;
                if sum(temp_adcs_pre)>0 && sum(temp_adcs_post)>0
                    ind1       = find(temp_adcs_pre==1);
                    ind1       = ind1(end); % last adc before trig
                    ind2       = find(temp_adcs_post==1);
                    ind2       = ind2(1);   % first adc after trig
                    k1         = sum(temp_adcs(1:ind1));
                    k2         = sum(temp_adcs(1:ind2));
                    dt_nominal = sum(SEQ.BLOCKS(ind1:ind2-1,2)) * dt;
                    dt_true    = adc_time_stamps(k2) - adc_time_stamps(k1);
                    dt_trig    = dt_true - dt_nominal;
                    SEQ.BLOCKS(ind_trig(j),2) = round(dt_trig/dt) + SEQ.BLOCKS(ind_trig(j),2);
                    SEQ.trig_delay(j,1)       = dt_trig;
                    SEQ.trig_timings(j,1)     = sum(SEQ.BLOCKS(1:ind_trig(j),2))*dt;                    
                else
                    SEQ.trig_delay(j,1)   = NaN;
                    SEQ.trig_timings(j,1) = NaN; 
                end
                clear temp_adcs temp_adcs_pre temp_adcs_post ind1 ind2 k1 k2 dt_nominal dt_true;
            end
            if isnan(SEQ.trig_timings(1)) && n_trig>1
                dt_nominal = sum(SEQ.BLOCKS(ind_trig(1):ind_trig(2),2)) * dt;
                SEQ.trig_timings(1) = SEQ.trig_timings(2) - dt_nominal;
            end
        end        
    end
end

%% -------------------- delete unnecessary blocks in case of kz partitions --------------------
function SEQ = filter_kz_partitions(SEQ, kz_part)
    if ~isempty(kz_part)
        if kz_part == 1
            for j = 1:size(SEQ.BLOCKS,1)
                if SEQ.BLOCKS(j,8)>0
                    temp_ext_id   = SEQ.BLOCKS(j, 8);
                    temp_ext_type = SEQ.EXTENSIONS(temp_ext_id, 2);
                    temp_ext_no   = SEQ.EXTENSIONS(temp_ext_id, 3);
                    if strcmp(SEQ.EXT_SPECS(temp_ext_type).type, 'LABELSET')
                        temp_ext_vals = strsplit(SEQ.EXT_SPECS(temp_ext_type).vals(temp_ext_no), ' ');
                        if strcmp(temp_ext_vals(3), 'START')
                            ind_start = j;
                        end
                        if strcmp(temp_ext_vals(3), 'STOP')
                            ind_stop = j;
                        end
                    end
                end
            end
            SEQ.BLOCKS = SEQ.BLOCKS(ind_start+1:ind_stop-1,:);
            SEQ.BLOCKS(:,1) = 1:size(SEQ.BLOCKS,1);
        end
    end
end

%% -------------------- correctt timgins for GE sequences --------------------
function SEQ = correct_GE_timings(SEQ)
    for j = 1:size(SEQ.BLOCKS,1)
        if SEQ.BLOCKS(j,8)>0
            temp_ext_id   = SEQ.BLOCKS(j, 8);
            temp_ext_type = SEQ.EXTENSIONS(temp_ext_id, 2);
            temp_ext_no   = SEQ.EXTENSIONS(temp_ext_id, 3);
            if strcmp(SEQ.EXT_SPECS(temp_ext_type).type, 'LABELSET')
                temp_ext_vals = strsplit(SEQ.EXT_SPECS(temp_ext_type).vals(temp_ext_no), ' ');
                if strcmp(temp_ext_vals(3), 'TRID')
                    if SEQ.BLOCKS(j,2) ~= 0
                        error('TRID label combined with delay!');
                    end
                    SEQ.BLOCKS(j,2) = round(117*1e-6 / SEQ.DEFINITIONS.blockRasterTime); % add approximately 117us delay for each TRID segment
                    SEQ.BLOCKS(j,8) = 0;
                end
            end
        end
    end
end

%% -------------------- delete unnecessary blocks in case of noise pre-scans ----------
function SEQ = delete_unnecessary_blocks(SEQ)

    % find input triggers
    temp_trig = zeros(size(SEQ.BLOCKS,1),1);
    for j = 1:size(SEQ.BLOCKS,1)
        if SEQ.BLOCKS(j,8)>0
            temp_ext_id   = SEQ.BLOCKS(j, 8);
            temp_ext_type = SEQ.EXTENSIONS(temp_ext_id, 2);
            temp_ext_no   = SEQ.EXTENSIONS(temp_ext_id, 3);
            if strcmp(SEQ.EXT_SPECS(temp_ext_type).type, 'TRIGGERS')
                temp_ext_vals = SEQ.EXT_SPECS(temp_ext_type).vals(temp_ext_no);
                temp_ext_vals = str2num(temp_ext_vals);
                if temp_ext_vals(2) == 2
                    temp_trig(j) = 1;
                end
            end
            clear temp_ext_id temp_ext_type temp_ext_no temp_ext_vals;
        end
    end

    % delete blocks before first rf or first input trigger
    ind_first_rf = find( (SEQ.BLOCKS(:,3)+temp_trig) > 0 );
    ind_first_rf = ind_first_rf(1);
    if ind_first_rf>1
        dt_before = sum(SEQ.BLOCKS(1:ind_first_rf-1,2));
        if ind_first_rf>2
            SEQ.BLOCKS(2:ind_first_rf-1,:) = [];
        end
        SEQ.BLOCKS(1,:) = [0 dt_before 0 0 0 0 0 0]; % conserve sequence duration
    end

    % delete blocks after last adc
    ind_last_adc = find(SEQ.BLOCKS(:,7)>0);
    ind_last_adc = ind_last_adc(end);
    if size(SEQ.BLOCKS,1)>ind_last_adc
        dt_after = sum(SEQ.BLOCKS(ind_last_adc+1:end,2));
        if size(SEQ.BLOCKS,1)-ind_last_adc > 1
            SEQ.BLOCKS(ind_last_adc+1:end-1,:) = [];
        end
        SEQ.BLOCKS(end,:) = [0 dt_after 0 0 0 0 0 0]; % conserve sequence duration
    end

    % repair block numbers
    SEQ.BLOCKS(:,1) = 1:size(SEQ.BLOCKS,1);

end

%% -------------------- calculate waveforms of sequence blocks --------------------
function SEQ = calc_block_waveforms(SEQ)

N = size(SEQ.BLOCKS,1);
cell_dt         = cell(N, 1);
cell_rf         = cell(N, 1);
cell_rf_use     = cell(N, 1);
cell_gx         = cell(N, 1);
cell_gy         = cell(N, 1);
cell_gz         = cell(N, 1);
cell_adc_on_off = cell(N, 1);
cell_adc_phase  = cell(N, 1);

dt_block = SEQ.DEFINITIONS.blockRasterTime;
dt_grad  = SEQ.DEFINITIONS.gradRasterTime;
dt_rf    = SEQ.DEFINITIONS.rfRasterTime;

parfor j = 1 : N
    temp_block = SEQ.BLOCKS(j,:);
    temp_dt    = temp_block(2) * dt_block;
    n          = round(temp_dt / 1e-6);

    %% calc rf waveform
    temp_rf_use = [];
    if temp_block(3) == 0
        temp_rf = zeros(n, 1);
    else
        % <id> <amp> <mag id> <phase id> <time id> <center> <delay> <freq_ppm> <phase_ppm> <freq_off> <phase_off> <use>
        %      <Hz>                                <us>     <us>    <ppm>      <rad/MHz>   <Hz>       <rad>
        temp_rf_specs = SEQ.RF(find(SEQ.RF(:,1)==temp_block(3)),:);
        temp_amp      = SEQ.SHAPES(temp_rf_specs(3)).data;
        temp_phase    = SEQ.SHAPES(temp_rf_specs(4)).data;
        if numel(temp_amp)==2 && isscalar(unique(temp_amp)) % block-pulse or spin-lock case
            temp_n   = SEQ.SHAPES(temp_rf_specs(5)).data;
            temp_n   = temp_n(2);
            temp_amp = ones(temp_n, 1) * temp_amp(1);
        end
        if numel(temp_phase)==2 && isscalar(unique(temp_phase)) % block-pulse or spin-lock case
            temp_n     = SEQ.SHAPES(temp_rf_specs(5)).data;
            temp_n     = temp_n(2);
            temp_phase = ones(temp_n, 1) * temp_phase(1);
        end        
        temp_amp     = interp1( linspace(0, 1, numel(temp_amp))',   temp_amp,   linspace(0, 1, round(numel(temp_amp)*dt_rf/1e-6))' );
        temp_phase   = interp1( linspace(0, 1, numel(temp_phase))', temp_phase, linspace(0, 1, round(numel(temp_phase)*dt_rf/1e-6))' );
        temp_t       = (1:numel(temp_amp))' * 1e-6;
        temp_rf      = temp_rf_specs(2) .* ... % amplitude [Hz]
                       temp_amp .* ... % shape [0...1]
                       exp(1i*2*pi*temp_phase) .* ... % phase [0...1]
                       exp(1i*2*pi*(SEQ.f0 * temp_rf_specs(8) * 1e6 + temp_rf_specs(10)) * temp_t) .* ... % offresonance [Hz]                  
                       exp(1i*(SEQ.f0 * temp_rf_specs(9) * 1e6 + temp_rf_specs(11))); % phase offset [0...2pi]                  
        temp_rf     = [zeros(temp_rf_specs(7),1); temp_rf];
        temp_rf     = [temp_rf; zeros(n-numel(temp_rf),1)];
        temp_rf_use = temp_rf_specs(12);
    end

    %% calc gx gy gz waveforms
    temp_gx = [];
    temp_gy = [];
    temp_gz = [];
    for xyz = 1 : 3
        if temp_block(3+xyz) == 0
            temp_g = zeros(n, 1);
        else
            temp_grad_specs = SEQ.GRADIENTS(find(SEQ.GRADIENTS(:,1)==temp_block(3+xyz)),:);
            if isempty(temp_grad_specs)
                temp_grad_specs = SEQ.TRAP(find(SEQ.TRAP(:,1)==temp_block(3+xyz)),:);
            end
            if numel(temp_grad_specs)==6
                % case: trapezoid gradient
                % <id> <amp>  <rise> <flat> <fall> <delay>
                %      <Hz/m> <us>   <us>   <us>   <us>
                temp_g = temp_grad_specs(2) * [zeros(temp_grad_specs(6),1); linspace(0,1,temp_grad_specs(3))'; ones(temp_grad_specs(4),1); linspace(1,0,temp_grad_specs(5))'];
                temp_g = [temp_g; zeros(n-numel(temp_g),1)]; 
            else
                % case: arbitrary gradient
                % <id> <amp>  <first> <last> <shape id> <time id> <delay>
                %      <Hz/m> <Hz/m>  <Hz/m>                      <us>
                temp_amp = SEQ.SHAPES(temp_grad_specs(5)).data;
                temp_amp = interp1( linspace(0, 1, numel(temp_amp))',   temp_amp,   linspace(0, 1, numel(temp_amp)*round(dt_grad/1e-6))' );
                temp_g   = temp_grad_specs(2) * temp_amp;
                temp_g   = [zeros(temp_grad_specs(7),1); temp_g];
                temp_g   = [temp_g; zeros(n-numel(temp_g),1)];
            end            
        end
        switch xyz
            case 1
                temp_gx = temp_g;
            case 2
                temp_gy = temp_g;
            case 3
                temp_gz = temp_g;
        end
    end
    % v1.5.1 -> here we should search for a rotation extension and apply
    % the quaternions to the waveforms... they are stored in SEQ.EXT_SPECS

    %% find adc sampling points
    if temp_block(7) == 0
        temp_adc_on_off = zeros(n, 1);
        temp_adc_phase  = zeros(n, 1);
    else
        % <id> <num> <dwell> <delay> <freq_ppm> <phase_ppm> <freq_off> <phase_off> <phase_shpe_id>
        %            <ns>    <us>    <ppm>      <rad/MHz>   <Hz>       <rad>
        temp_adc_specs  = SEQ.ADC(find(SEQ.ADC(:,1)==temp_block(7)),:);
        temp_adc_on_off = [zeros(temp_adc_specs(4),1); ones(round(temp_adc_specs(2) * temp_adc_specs(3) *1e-9 /1e-6),1)];        
        temp_adc_on_off = [temp_adc_on_off; zeros(n-numel(temp_adc_on_off),1)];
        temp_adc_phase  = ones(n, 1) * (SEQ.f0 * temp_adc_specs(6) *1e6 + temp_adc_specs(8));
    end

    %% interpolate to target raster
    if SEQ.dt ~= 1e-6
        temp_rf         = interp1(linspace(0, 1, n)', temp_rf,         linspace(0, 1, round(n*1e-6/SEQ.dt))');
        temp_gx         = interp1(linspace(0, 1, n)', temp_gx,         linspace(0, 1, round(n*1e-6/SEQ.dt))');
        temp_gy         = interp1(linspace(0, 1, n)', temp_gy,         linspace(0, 1, round(n*1e-6/SEQ.dt))');
        temp_gz         = interp1(linspace(0, 1, n)', temp_gz,         linspace(0, 1, round(n*1e-6/SEQ.dt))');
        temp_adc_on_off = interp1(linspace(0, 1, n)', temp_adc_on_off, linspace(0, 1, round(n*1e-6/SEQ.dt))');
        temp_adc_phase  = interp1(linspace(0, 1, n)', temp_adc_phase,  linspace(0, 1, round(n*1e-6/SEQ.dt))');
    end

    %% save to cell array
    cell_dt{j}         = temp_dt;
    cell_rf{j}         = temp_rf;
    cell_rf_use{j}     = temp_rf_use;
    cell_gx{j}         = temp_gx;
    cell_gy{j}         = temp_gy;
    cell_gz{j}         = temp_gz;
    cell_adc_on_off{j} = temp_adc_on_off;
    cell_adc_phase{j}  = wrapTo2Pi(temp_adc_phase);

end

% save to global struct
SEQ.WAVES = struct();
for j = 1:N
    SEQ.WAVES(j,1).dt         = cell_dt{j};
    SEQ.WAVES(j,1).RF         = cell_rf{j};
    SEQ.WAVES(j,1).RF_use     = cell_rf_use{j};
    SEQ.WAVES(j,1).GX         = cell_gx{j};
    SEQ.WAVES(j,1).GY         = cell_gy{j};
    SEQ.WAVES(j,1).GZ         = cell_gz{j};
    SEQ.WAVES(j,1).ADC_ON_OFF = cell_adc_on_off{j};
    SEQ.WAVES(j,1).ADC_PHASE  = cell_adc_phase{j};
end

end

%% -------------------- find echo positions --------------------
function SEQ = find_echo_positions(SEQ, echo_mode, adc_npad)

    % this functions searches for echo positions within the ADCs
    % select echo mode:
    % 'spiral_out' -> the first ADC point is used as the echo position
    % 'center'     -> the center of the ADC is used as the echo position
    % 'auto'       -> calculate kx, ky, kz and search for local minima (use e.g. for rosettes)
    
    if strcmp(echo_mode, 'auto')
        SEQ = calc_kxyz(SEQ);
    end
    
    for j = 1:size(SEQ.WAVES, 1)
        if sum(SEQ.WAVES(j).ADC_ON_OFF)==0
            SEQ.WAVES(j).ADC_SIM = zeros(size(SEQ.WAVES(j).ADC_ON_OFF));
        else
            switch echo_mode
                case 'spiral_out'
                    temp_adc   = SEQ.WAVES(j).ADC_ON_OFF;
                    temp_start = find(temp_adc==1);
                    temp_start = temp_start(1);
                    if ~isempty(adc_npad)
                        temp_start = temp_start + adc_npad(1);
                    end                    
                    temp_adc   = temp_adc*0;
                    temp_adc(temp_start) = 1;
                    SEQ.WAVES(j).ADC_SIM = temp_adc;
                    clear temp_adc temp_start;
                case 'center'
                    temp_adc   = SEQ.WAVES(j).ADC_ON_OFF;
                    temp_start = find(temp_adc==1);
                    temp_start = temp_start(1);
                    temp_end   = find(temp_adc==1);
                    temp_end   = temp_end(end);
                    temp_adc   = temp_adc*0;
                    temp_adc(round((temp_start+temp_end)/2)) = 1;
                    SEQ.WAVES(j).ADC_SIM = temp_adc;
                    clear temp_adc temp_start temp_end;
                case 'auto'
                    temp_kxyz = SEQ.WAVES(j).KX.^2+SEQ.WAVES(j).KY.^2+SEQ.WAVES(j).KZ.^2;
                    temp_adc  = SEQ.WAVES(j).ADC_ON_OFF;
                    k = 1;
                    while temp_adc(k)==0
                        temp_kxyz(1) = [];
                        k = k + 1;
                    end
                    temp_adcpad = numel(temp_adc)- numel(temp_kxyz);
                    k = numel(temp_adc);
                    while temp_adc(k)==0
                        temp_kxyz(end) = [];
                        k = k - 1;
                    end
                    temp_kxyz = temp_kxyz - min(temp_kxyz);                
                    temp1      = temp_kxyz < 0.005*median(temp_kxyz);
                    temp2      = diff([0; temp1; 0]);
                    temp_start = find(temp2 == 1);
                    temp_end   = find(temp2 == -1) - 1;
                    temp_pos   = round((temp_start + temp_end) / 2);
                    if temp1(1)==1
                        temp_pos(1) = 1;
                    end
                    temp_pos = temp_pos + temp_adcpad;
                    temp_adc = temp_adc*0;
                    temp_adc(temp_pos) = 1;
                    SEQ.WAVES(j).ADC_SIM = temp_adc;
                    clear temp_kxyz temp_adc k temp_adcpad temp1 temp2 temp_start temp_end temp_pos;
            end
        end
    end
end

%% -------------------- calculate kx ky kz --------------------
function SEQ = calc_kxyz(SEQ)
    kx = 0;
    ky = 0;
    kz = 0;    
    for j = 1:size(SEQ.WAVES, 1)
    
        if isempty(SEQ.WAVES(j).RF_use)
            if isempty(SEQ.WAVES(j).ADC_ON_OFF)
                if ~isempty(SEQ.WAVES(j).GX)
                    kx = kx + sum(SEQ.WAVES(j).GX)*SEQ.dt;
                end
                if ~isempty(SEQ.WAVES(j).GY)
                    ky = ky + sum(SEQ.WAVES(j).GY)*SEQ.dt;
                end
                if ~isempty(SEQ.WAVES(j).GZ)
                    kz = kz + sum(SEQ.WAVES(j).GZ)*SEQ.dt;
                end            
            else
                if isempty(SEQ.WAVES(j).GX)
                    temp_kx = kx + zeros(numel(SEQ.WAVES(j).GX),1);
                else
                    temp_kx = kx + cumsum(SEQ.WAVES(j).GX)*SEQ.dt;
                    kx = kx + sum(SEQ.WAVES(j).GX)*SEQ.dt;                
                end
                if isempty(SEQ.WAVES(j).GY)
                    temp_ky = ky + zeros(numel(SEQ.WAVES(j).GY),1);
                else
                    temp_ky = ky + cumsum(SEQ.WAVES(j).GY)*SEQ.dt;
                    ky = ky + sum(SEQ.WAVES(j).GY)*SEQ.dt;                
                end
                if isempty(SEQ.WAVES(j).GZ)
                    temp_kz = kz + zeros(numel(SEQ.WAVES(j).GZ),1);
                else
                    temp_kz = kz + cumsum(SEQ.WAVES(j).GZ)*SEQ.dt;
                    kz = kz + sum(SEQ.WAVES(j).GZ)*SEQ.dt;                
                end            
                SEQ.WAVES(j).KX = temp_kx;
                SEQ.WAVES(j).KY = temp_ky;
                SEQ.WAVES(j).KZ = temp_kz;
                clear temp_kx temp_ky temp_kz;
            end
        else
            temp_ind = round(median(find(abs(SEQ.WAVES(j).RF)==max(abs(SEQ.WAVES(j).RF)))));
            if ~isempty(SEQ.WAVES(j).GX)
                kx = kx + sum(SEQ.WAVES(j).GX(1:temp_ind))*SEQ.dt;
            end
            if ~isempty(SEQ.WAVES(j).GY)
                ky = ky + sum(SEQ.WAVES(j).GY(1:temp_ind))*SEQ.dt;
            end
            if ~isempty(SEQ.WAVES(j).GZ)
                kz = kz + sum(SEQ.WAVES(j).GZ(1:temp_ind))*SEQ.dt;
            end
            switch SEQ.WAVES(j).RF_use
                case 2 % refocusing
                    kx = -kx;
                    ky = -ky;
                    kz = -kz;
                otherwise % excitation
                    kx = 0;
                    ky = 0;
                    kz = 0;
            end
            if ~isempty(SEQ.WAVES(j).GX)
                kx = kx + sum(SEQ.WAVES(j).GX(temp_ind+1:end))*SEQ.dt;
            end
            if ~isempty(SEQ.WAVES(j).GY)
                ky = ky + sum(SEQ.WAVES(j).GY(temp_ind+1:end))*SEQ.dt;
            end
            if ~isempty(SEQ.WAVES(j).GZ)
                kz = kz + sum(SEQ.WAVES(j).GZ(temp_ind+1:end))*SEQ.dt;
            end        
        end    
    end
end

%% -------------------- calculate full sequence waveforms --------------------
function SEQ = calc_full_waveforms(SEQ, dt, flag_plot)

    % calculate array size
    t_end = sum(SEQ.BLOCKS(:,2)) * SEQ.DEFINITIONS.blockRasterTime;
    N     = round(t_end/dt);
    
    % init waveform arrays
    RF         = double(zeros(N, 1)); % [Hz]
    GX         = double(zeros(N, 1)); % [Hz/m]
    GY         = double(zeros(N, 1)); % [Hz/m]
    GZ         = double(zeros(N, 1)); % [Hz/m]
    ADC_ON_OFF = double(zeros(N, 1)); % [0, 1]
    ADC_PHASE  = double(zeros(N, 1)); % [rad]
    ADC_SIM    = double(zeros(N, 1)); % [0 1]
    
    ind_start = 1;
    for j = 1:size(SEQ.BLOCKS,1)
        n = numel(SEQ.WAVES(j).RF);
        RF(ind_start:ind_start+n-1)         = SEQ.WAVES(j).RF;
        GX(ind_start:ind_start+n-1)         = SEQ.WAVES(j).GX;
        GY(ind_start:ind_start+n-1)         = SEQ.WAVES(j).GY;
        GZ(ind_start:ind_start+n-1)         = SEQ.WAVES(j).GZ;
        ADC_ON_OFF(ind_start:ind_start+n-1) = SEQ.WAVES(j).ADC_ON_OFF;
        ADC_PHASE(ind_start:ind_start+n-1)  = SEQ.WAVES(j).ADC_PHASE;
        ADC_SIM(ind_start:ind_start+n-1)    = SEQ.WAVES(j).ADC_SIM;
        ind_start = ind_start + n;
    end
    clear ind_start j n;

    % plot waveforms
    if flag_plot==2
        t = (1:N)'*dt*1e3; % [ms]
        
        fig = figure();
        ax1 = subplot(3,2,1);
        hold on
        xline(find(ADC_SIM==1)*dt*1e3, 'k-');
        plot(t(ADC_ON_OFF==1), ADC_ON_OFF(ADC_ON_OFF==1), 'rx');
        xlim([t(1) t(end)]);
        ylabel('ADC on/off');
        yticks(1);
        set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
        
        ax3 = subplot(3,2,3);
        plot(t, abs(RF), 'b-');
        xlim([t(1) t(end)]);
        ylabel('RF magnitude [Hz]');
        set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
        
        ax5 = subplot(3,2,5);
        hold on
        plot(t, angle(RF), 'b-');
        plot(t(ADC_ON_OFF==1), ADC_PHASE(ADC_ON_OFF==1), 'rx');
        xlim([t(1) t(end)]);
        ylabel('RF/ADC phase [rad]');
        xlabel('time [ms]')
        set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
        
        ax2 = subplot(3,2,2);
        plot(t, GX/1e3, 'k-');
        xlim([t(1) t(end)]);
        ylabel('x gradient [kHz/m]');
        set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
        
        ax4 = subplot(3,2,4);
        plot(t, GY/1e3, 'k-');
        xlim([t(1) t(end)]);
        ylabel('y gradient [kHz/m]');
        set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
        
        ax6 = subplot(3,2,6);
        plot(t, GZ/1e3, 'k-');
        xlim([t(1) t(end)]);
        ylabel('z gradient [kHz/m]');
        set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
        
        linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')
        h = zoom(fig);
        setAxesZoomMotion(h, ax1, 'horizontal');
        
        if isfield(SEQ.DEFINITIONS, 'name')
            sgtitle(SEQ.DEFINITIONS.name, 'Interpreter', 'none', 'FontWeight', 'bold');
        end
    end

    % store arrays in SEQ struct
    SEQ.FULL.RF         = RF;
    SEQ.FULL.GX         = GX;
    SEQ.FULL.GY         = GY;
    SEQ.FULL.GZ         = GZ;
    SEQ.FULL.ADC_ON_OFF = ADC_ON_OFF;
    SEQ.FULL.ADC_PHASE  = ADC_PHASE;
    SEQ.FULL.ADC_SIM    = ADC_SIM;

end

%% -------------------- reduce memory of seq waves struct --------------------
function SEQ = reduce_seq_wave_size(SEQ)
    for j = 1:size(SEQ.BLOCKS,1)
        if SEQ.BLOCKS(j,3)==0
            SEQ.WAVES(j).RF = [];
        end
        if SEQ.BLOCKS(j,4)==0
            SEQ.WAVES(j).GX = [];
        end
        if SEQ.BLOCKS(j,5)==0
            SEQ.WAVES(j).GY = [];
        end
        if SEQ.BLOCKS(j,6)==0
            SEQ.WAVES(j).GZ = [];
        end
        if SEQ.BLOCKS(j,7)==0
            SEQ.WAVES(j).ADC_ON_OFF = [];
            SEQ.WAVES(j).ADC_PHASE  = [];
            SEQ.WAVES(j).ADC_SIM    = [];
        end
    end
end

%% -------------------- init operator list for brute force bloch simulation --------------------
function SIM = init_brute_force_bloch_simulation(SEQ)

% ----- versions: -----
% V1: 17.04.2025; M. Gram; University of Wuerzburg
% V2: 15.06.2025: M. Gram; University of Wuerzburg; adiabatic spin-locking

% this function generates operator lists for the simulation of .seq files
% the operator lists and corresponding waveforms (RF, GX, GY, ...) are
% compressed for eliminating unnecessary simulation steps
% e.g. multiple delays or gradient moments are joined
% the operator list can be simulated with MRF_sim_brute_force()

% ----- input: -----
% SEQ.WAVES:  struct, which contains the full .seq information
%   .RF:           [Hz]   complex valued waveform of the radiofrequency pulses
%   .GX, .GY, .GZ: [Hz/m] waveform of the x,y,z gradients
%   .ADC_SIM:      [0 1]  time stamps of simulated acquisitions
%   .ADC_PHASE:    [rad]  relative phase of receiver

% ----- output: -----
% SIM
%   .ID:  compressed ID array for selecting the correct Bloch operator
%   .RF:  compressed complex RF waveforms [rad/s]
%   .GX:  compressed x gradient moments [rad/m] or waveforms [rad/s/m]
%   .GY:  compressed y gradient moments [rad/m] or waveforms [rad/s/m]
%   .GZ:  compressed z gradient moments [rad/m] or waveforms [rad/s/m]
%   .DB:  compressed b-value of the diffusion encoding tensor [s/mm^2]
%   .DT:  compressed variable time steps [s]
%   .PHI: adc reciever phases [rad]

% ----- ID: -----
% 0:  adc
% 1:  free relaxation & dephasing
% 2:  rf: global standard (e.g. inversion, saturation, fat suppression)
% 3:  rf: global adiabatic (e.g. AHP, HYPSECH)
% 4:  rf: global spin-lock (continuous wave block pulse)
% 5:  rf: slice selective excitation (e.g. sinc, gauss, slr)
% 6:  rf: slice selective refocusing (e.g. sinc, gauss, slr) !!! to do
% 7:  rf: slice selective adiabatic (e.g. goia-wurst) !!! to do
% 8:  gradient dephasing: x
% 9:  gradient dephasing: y
% 10: gradient dephasing: z
% 11: diffusion encoding
% 12: adiabatic spin-lock

%% approx maximum size of operator lists
max_idx = 0;
for j=1:numel(SEQ.WAVES)
    if ~isempty(SEQ.WAVES(j).RF)
        max_idx = max_idx + numel(SEQ.WAVES(j).RF) + 10;
    else
        max_idx = max_idx + 10;
    end
end

% init operator lists
ID = zeros(max_idx, 1) - 1;
RF = zeros(max_idx, 1);
GX = zeros(max_idx, 1);
GY = zeros(max_idx, 1);
GZ = zeros(max_idx, 1);
DB = zeros(max_idx, 1);
DT = zeros(max_idx, 1);
clear max_idx;

%% sequence blocks -> ID list and waveforms

SEQ.BLOCKS = SEQ.BLOCKS>0;
SEQ.BLOCKS(:,1:2) = [];
% rf gx gy gz adc (ext)
% [SEQ_BLOCKS ] -> operator
% [0 0 0 0 0 x] -> free or trigger
% [1 0 0 1 0 x] -> slice selective rf, gz
% [1 0 0 0 0 x] -> global rf, adiabatic, fat, spin-lock
% [0 x y z 1 x] -> gx, gy, gz, adc
% [0 x y z 0 x] -> gx gy gz
% note: the ext block can be ignored, since trigger timings were already
% corrected and label statements are not important for the simulation

j_adc = 1;
idx   = 1;

%% scan all sequence blocks
for j = 1:size(SEQ.BLOCKS,1)

    temp_block = SEQ.BLOCKS(j,1:5);

    %% adc
    if temp_block(5)==1

        % read gradients and adcs
        temp_gx         = SEQ.WAVES(j).GX;
        temp_gy         = SEQ.WAVES(j).GY;
        temp_gz         = SEQ.WAVES(j).GZ;
        temp_adc_on_off = SEQ.WAVES(j).ADC_SIM;
        temp_adc_on_off = find(temp_adc_on_off==1);

        % free relaxation before first adc
        temp_range = 1 : temp_adc_on_off(1)-1;
        if ~isempty(temp_range)
            ID(idx) = 1;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = numel(temp_range) *SEQ.dt;
            idx = idx + 1;
        end
        % x dephasing before first adc
        if ~isempty(temp_gx)
        if abs(sum(temp_gx(temp_range))*SEQ.dt) > 0
            ID(idx) = 8;
            RF(idx) = 0;
            GX(idx) = sum(temp_gx(temp_range))*SEQ.dt;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end
        end
        % y dephasing before first adc
        if ~isempty(temp_gy)
        if abs(sum(temp_gy(temp_range))*SEQ.dt) > 0
            ID(idx) = 9;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = sum(temp_gy(temp_range))*SEQ.dt;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end
        end
        % z dephasing before first adc
        if ~isempty(temp_gz)
        if abs(sum(temp_gz(temp_range))*SEQ.dt) > 0
            ID(idx) = 10;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = sum(temp_gz(temp_range)) *SEQ.dt;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end
        end
        % diffusion encoding before first adc
        ID(idx) = 11;
        RF(idx) = 0;
        GX(idx) = 0;
        GY(idx) = 0;
        GZ(idx) = 0;
        DB(idx) = calc_diffusion_b_value(SEQ.dt, temp_gx, temp_gy, temp_gz, temp_range);
        DT(idx) = 0;
        idx = idx + 1;

        % ADCs and free relaxation between ADCs
        temp_adc_on_off = [temp_adc_on_off; round(SEQ.WAVES(j).dt/SEQ.dt)];
        for k = 1:numel(temp_adc_on_off)-1    
            temp_range = temp_adc_on_off(k) : temp_adc_on_off(k+1);
            % adc
            ID(idx) = 0;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;            
            PHI(j_adc, 1) = SEQ.WAVES(j).ADC_PHASE(temp_adc_on_off(k));
            j_adc = j_adc + 1;
            idx = idx + 1;
            % free relaxation between ADCs
            ID(idx) = 1;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = numel(temp_range) *SEQ.dt;
            idx = idx + 1;
        end        

        % remaining x dephasing
        temp_range = temp_adc_on_off(1) : temp_adc_on_off(end);
        if ~isempty(temp_gx)
        if abs(sum(temp_gx(temp_range))*SEQ.dt) > 0
            ID(idx) = 8;
            RF(idx) = 0;
            GX(idx) = sum(temp_gx(temp_range))*SEQ.dt;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end
        end
        % remaining y dephasing
        if ~isempty(temp_gy)
        if abs(sum(temp_gy(temp_range))*SEQ.dt) > 0
            ID(idx) = 9;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = sum(temp_gy(temp_range))*SEQ.dt;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end
        end
        % remaining z dephasing
        if ~isempty(temp_gz)
        if abs(sum(temp_gz(temp_range))*SEQ.dt) > 0
            ID(idx) = 10;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = sum(temp_gz(temp_range)) *SEQ.dt;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end
        end
        % remaining diffusion encoding
        ID(idx) = 11;
        RF(idx) = 0;
        GX(idx) = 0;
        GY(idx) = 0;
        GZ(idx) = 0;
        DB(idx) = calc_diffusion_b_value(SEQ.dt, temp_gx, temp_gy, temp_gz, temp_range);
        DT(idx) = 0;
        idx = idx + 1;
        clear temp_adc_on_off temp_gx temp_gy temp_gz temp_range;

    end
    
    %% free relaxation
    if all(temp_block == [0 0 0 0 0])
        ID(idx) = 1;
        RF(idx) = 0;
        GX(idx) = 0;
        GY(idx) = 0;
        GZ(idx) = 0;
        DB(idx) = 0;
        DT(idx) = SEQ.WAVES(j).dt;
        idx = idx + 1;
    end    

    %% global excitation
    if all(temp_block == [1 0 0 0 0])
        temp_rf = SEQ.WAVES(j).RF;

        % free relaxation during rf delay
        while temp_rf(1)==0
            temp_rf(1) = [];
        end
        if numel(temp_rf) < numel(SEQ.WAVES(j).RF)
            ID(idx) = 1;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = (numel(SEQ.WAVES(j).RF) - numel(temp_rf)) *SEQ.dt;
            idx = idx + 1;           
        end
        temp_n_ringdown = numel(temp_rf);
        while temp_rf(end)==0
            temp_rf(end) = [];
        end
        temp_n_ringdown = temp_n_ringdown - numel(temp_rf);

        % test for adiabatic pulses
        temp_off = abs(diff(unwrap(angle(temp_rf))));
        temp_off = unique(temp_off);

        % adiabatic spin-lock
        if SEQ.WAVES(j).RF_use==6
            idx     = idx + (1:numel(temp_rf)) - 1;
            ID(idx) = ones(numel(temp_rf),1)*12;
            RF(idx) = temp_rf;
            GX(idx) = zeros(numel(temp_rf),1);
            GY(idx) = zeros(numel(temp_rf),1);
            GZ(idx) = zeros(numel(temp_rf),1);
            DB(idx) = zeros(numel(temp_rf),1);
            DT(idx) = ones(numel(temp_rf),1)*SEQ.dt;
            idx     = idx(end) + 1;

        % continuous wave -> spin-locking
        elseif isscalar(unique(abs(temp_rf)))
            temp_switch = find(diff((angle(temp_rf)))~=0);
            temp_switch = [0; temp_switch; numel(temp_rf)];
            for k=1:numel(temp_switch)-1
                temp_sl = temp_rf(temp_switch(k)+1:temp_switch(k+1));
                ID(idx) = 4;
                RF(idx) = temp_sl(1);
                GX(idx) = 0;
                GY(idx) = 0;
                GZ(idx) = 0;
                DB(idx) = 0;
                DT(idx) = numel(temp_sl)*SEQ.dt;
                idx = idx + 1;
            end
            clear k temp_sl temp_switch;

        % adiabatic pulse
        elseif numel(temp_off) > round(numel(temp_rf)*0.05) % more than 5% phase modulation
            idx     = idx + (1:numel(temp_rf)) - 1;
            ID(idx) = ones(numel(temp_rf),1)*3;
            RF(idx) = temp_rf;
            GX(idx) = zeros(numel(temp_rf),1);
            GY(idx) = zeros(numel(temp_rf),1);
            GZ(idx) = zeros(numel(temp_rf),1);
            DB(idx) = zeros(numel(temp_rf),1);
            DT(idx) = ones(numel(temp_rf),1)*SEQ.dt;
            idx     = idx(end) + 1;

        % other global excitation pulse: inversion, saturation, fat suppression
        else
            idx     = idx + (1:numel(temp_rf)) - 1;
            ID(idx) = ones(numel(temp_rf),1)*2;
            RF(idx) = temp_rf;
            GX(idx) = zeros(numel(temp_rf),1);
            GY(idx) = zeros(numel(temp_rf),1);
            GZ(idx) = zeros(numel(temp_rf),1);
            DB(idx) = zeros(numel(temp_rf),1);
            DT(idx) = ones(numel(temp_rf),1)*SEQ.dt;
            idx     = idx(end) + 1;
        end

        % free relaxation during rf ringdown
        if temp_n_ringdown>0
            ID(idx) = 1;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = temp_n_ringdown*SEQ.dt;
            idx = idx + 1;            
        end

        clear temp_n_ringdown temp_rf temp_off;

    end

    %% slice selective excitation
    if all(temp_block == [1 0 0 1 0])
        temp_rf = SEQ.WAVES(j).RF;
        temp_gz = SEQ.WAVES(j).GZ;

        % free relaxation and z dephasing and diffusion encoding before start of rf pulse
        while temp_rf(1)==0
            temp_rf(1) = [];
        end
        if numel(temp_gz) > numel(temp_rf)
            % free relaxation
            ID(idx) = 1;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = (numel(temp_gz)-numel(temp_rf))*SEQ.dt;
            idx = idx + 1;
            % z dephasing
            ID(idx) = 10;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = sum(temp_gz(1:(numel(temp_gz)-numel(temp_rf))))*SEQ.dt;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
            % diffusion encoding
            ID(idx) = 11;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = calc_diffusion_b_value(SEQ.dt, [], [], temp_gz(1:(numel(temp_gz)-numel(temp_rf))));
            DT(idx) = 0;
            idx = idx + 1;          
            temp_gz(1:(numel(temp_gz)-numel(temp_rf))) = [];
        end

        % rf pulse
        while temp_rf(end)==0
            temp_rf(end) = [];
        end
        temp_off = abs(diff(unwrap(angle(temp_rf))));
        temp_off = unique(temp_off);
        idx      = idx + (1:numel(temp_rf)) - 1;
        if numel(temp_off) > round(numel(temp_rf)*0.05) % more than 5% phase modulation
            ID(idx) = ones(numel(temp_rf),1)*7; % slice selective adiabatic
        elseif (abs(sum(temp_rf))*SEQ.dt*360) > 135
            ID(idx) = ones(numel(temp_rf),1)*6; % slice selective refocusing
        else
            ID(idx) = ones(numel(temp_rf),1)*5; % slice selective excitation
        end
        RF(idx) = temp_rf;
        GX(idx) = zeros(numel(temp_rf),1);
        GY(idx) = zeros(numel(temp_rf),1);
        GZ(idx) = temp_gz(1:numel(temp_rf));
        DB(idx) = zeros(numel(temp_rf),1);
        DT(idx) = ones(numel(temp_rf),1)*SEQ.dt;
        idx     = idx(end) + 1;

        % free relaxation and z dephasing and diffusion encoding after end of rf pulse
        if numel(temp_gz) > numel(temp_rf)
            % free relaxation
            ID(idx) = 1;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = (numel(temp_gz)-numel(temp_rf))*SEQ.dt;
            idx = idx + 1;
            % z dephasing
            ID(idx) = 10;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = sum(temp_gz(numel(temp_rf)+1:end))*SEQ.dt;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
            % diffusion encoding
            ID(idx) = 11;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = calc_diffusion_b_value(SEQ.dt, [], [], temp_gz(numel(temp_rf)+1:end));
            DT(idx) = 0;
            idx = idx + 1;            
        end

        clear temp_rf temp_gz temp_off;

    end

    %% only gradients
    if all([temp_block(1) temp_block(5)] == [0 0])
    if sum([temp_block(2) temp_block(3) temp_block(4)])>0
        temp_gx = SEQ.WAVES(j).GX;
        temp_gy = SEQ.WAVES(j).GY;
        temp_gz = SEQ.WAVES(j).GZ;

        % free relaxation
        ID(idx) = 1;
        RF(idx) = 0;
        GX(idx) = 0;
        GY(idx) = 0;
        GZ(idx) = 0;
        DB(idx) = 0;
        DT(idx) = SEQ.WAVES(j).dt;
        idx = idx + 1; 

        % x dephasing
        if abs(sum(temp_gx)*SEQ.dt) > 0
            ID(idx) = 8;
            RF(idx) = 0;
            GX(idx) = sum(temp_gx)*SEQ.dt;
            GY(idx) = 0;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end

        % y dephasing
        if abs(sum(temp_gy)*SEQ.dt) > 0
            ID(idx) = 9;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = sum(temp_gy)*SEQ.dt;
            GZ(idx) = 0;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end

        % z dephasing
        if abs(sum(temp_gz)*SEQ.dt) > 0
            ID(idx) = 10;
            RF(idx) = 0;
            GX(idx) = 0;
            GY(idx) = 0;
            GZ(idx) = sum(temp_gz) *SEQ.dt;
            DB(idx) = 0;
            DT(idx) = 0;
            idx = idx + 1;
        end       
        
        % diffusion encoding
        ID(idx) = 11;
        RF(idx) = 0;
        GX(idx) = 0;
        GY(idx) = 0;
        GZ(idx) = 0;
        DB(idx) = calc_diffusion_b_value(SEQ.dt, temp_gx, temp_gy, temp_gz);
        DT(idx) = 0;
        idx = idx + 1;
        clear temp_gx temp_gy temp_gz;

    end
    end

    clear temp_block;

end
clear j c;

% remove unnecessary entries
RF(ID==-1) = [];
GX(ID==-1) = [];
GY(ID==-1) = [];
GZ(ID==-1) = [];
DB(ID==-1) = [];
DT(ID==-1) = [];
ID(ID==-1) = [];

%% join relaxation, spoiler, crusher and diffusion weighting

% init joined arrays
ID_join     = zeros(size(ID)) - 1;
RF_join     = zeros(size(RF));
GX_join     = zeros(size(GX));
GY_join     = zeros(size(GY));
GZ_join     = zeros(size(GZ));
DB_join     = zeros(size(DB));
DT_join     = zeros(size(DT));
ID_join(1)  = ID(1);
RF_join(1)  = RF(1);
GX_join(1)  = GX(1);
GY_join(1)  = GY(1);
GZ_join(1)  = GZ(1);
DB_join(1)  = DB(1);
DT_join(1)  = DT(1);
j_start     = 1;
idx         = 2;
temp_rf_adc = ismember(ID, [0 2 3 4 5 6 7 12]);

while j_start<numel(ID)
    j_end = j_start + 1;
    while temp_rf_adc(j_end)==0
        j_end = j_end + 1;
        if j_end>numel(ID)
            j_end = numel(ID);
            break;
        end
    end
    if j_end > j_start+1
        temp_ind = j_start+1 : j_end-1;
        if abs(sum(DT(temp_ind))) > 0
            ID_join(idx) = 1;
            RF_join(idx) = 0;
            GX_join(idx) = 0;
            GY_join(idx) = 0;
            GZ_join(idx) = 0;
            DB_join(idx) = 0;
            DT_join(idx) = sum(DT(temp_ind));
            idx = idx + 1;
        end
        if abs(sum(GX(temp_ind))) > 0
            ID_join(idx) = 8;
            RF_join(idx) = 0;
            GX_join(idx) = sum(GX(temp_ind));
            GY_join(idx) = 0;
            GZ_join(idx) = 0;
            DB_join(idx) = 0;
            DT_join(idx) = 0;
            idx = idx + 1;
        end
        if abs(sum(GY(temp_ind))) > 0
            ID_join(idx) = 9;
            RF_join(idx) = 0;
            GX_join(idx) = 0;
            GY_join(idx) = sum(GY(temp_ind));
            GZ_join(idx) = 0;
            DB_join(idx) = 0;
            DT_join(idx) = 0;
            idx = idx + 1;
        end
        if abs(sum(GZ(temp_ind))) > 0
            ID_join(idx) = 10;
            RF_join(idx) = 0;
            GX_join(idx) = 0;
            GY_join(idx) = 0;
            GZ_join(idx) = sum(GZ(temp_ind));
            DB_join(idx) = 0;
            DT_join(idx) = 0;
            idx = idx + 1;
        end        
        if abs(sum(DB(temp_ind))) > 0
            ID_join(idx) = 11;
            RF_join(idx) = 0;
            GX_join(idx) = 0;
            GY_join(idx) = 0;
            GZ_join(idx) = 0;
            DB_join(idx) = sum(DB(temp_ind));
            DT_join(idx) = 0;
            idx = idx + 1;
        end
    end
    ID_join(idx) = ID(j_end);
    RF_join(idx) = RF(j_end);
    GX_join(idx) = GX(j_end);
    GY_join(idx) = GY(j_end);
    GZ_join(idx) = GZ(j_end);
    DB_join(idx) = DB(j_end);
    DT_join(idx) = DT(j_end);
    idx = idx + 1;
    j_start = j_end ;
end

% remove unnecessary entries
RF_join(ID_join==-1) = [];
GX_join(ID_join==-1) = [];
GY_join(ID_join==-1) = [];
GZ_join(ID_join==-1) = [];
DB_join(ID_join==-1) = [];
DT_join(ID_join==-1) = [];
ID_join(ID_join==-1) = [];

clear ID RF GX GY GZ DB DT j_start j_end temp_ind;

%% eliminate nearly zero x and y gradients; abs(moment) < 1 (this means: 1/m = 360°/m = 3.6°/cm) 
temp_del = (abs(GX_join)<1) .* (ID_join==8);
ID_join(temp_del==1) = [];
RF_join(temp_del==1) = [];
GX_join(temp_del==1) = [];
GY_join(temp_del==1) = [];
GZ_join(temp_del==1) = [];
DB_join(temp_del==1) = [];
DT_join(temp_del==1) = [];
temp_del = (abs(GY_join)<1) .* (ID_join==9);
ID_join(temp_del==1) = [];
RF_join(temp_del==1) = [];
GX_join(temp_del==1) = [];
GY_join(temp_del==1) = [];
GZ_join(temp_del==1) = [];
DB_join(temp_del==1) = [];
DT_join(temp_del==1) = [];
clear temp_del;

%% eliminate nearly zero diffusion weighting; if b<1 and ADC=0.005 -> exp(-b*ADC)=0.995
temp_del = (abs(DB_join)<1) .* (ID_join==11);
ID_join(temp_del==1) = [];
RF_join(temp_del==1) = [];
GX_join(temp_del==1) = [];
GY_join(temp_del==1) = [];
GZ_join(temp_del==1) = [];
DB_join(temp_del==1) = [];
DT_join(temp_del==1) = [];
clear temp_del;

%% output SIM struct, convert Hz -> rad/s
SIM.ID  = int32(ID_join);
SIM.RF  = RF_join *2*pi;    
SIM.GX  = GX_join *2*pi;
SIM.GY  = GY_join *2*pi;
SIM.GZ  = GZ_join *2*pi;
SIM.DB  = DB_join;
SIM.DT  = DT_join;
SIM.PHI = PHI;

end

%% --------------------- calculate diffusion b-value ---------------------
function b = calc_diffusion_b_value(dt, gx, gy, gz, idx)

% ----- reference: -----
% Gradient waveform design for tensor-valued encoding in diffusion MRI
% Filip Szczepankiewicz, Carl-Fredrik Westin, Markus Nilsson 
% Journal of Neuroscience Methods
% Volume 348, 15 January 2021, 109007
% doi.org/10.1016/j.jneumeth.2020.109007

% gx, gy, gz: column vectors [Hz/m]
% b: b-value of the diffusion encoding tensor [s/mm^2]

% ----- empty inputs and index range -----
if nargin==5
    if ~isempty(gx)
        gx = gx(idx);
    end
    if ~isempty(gy)
        gy = gy(idx);
    end
    if ~isempty(gz)
        gz = gz(idx);
    end
end
N = max([numel(gx); numel(gy); numel(gz)]);
if isempty(gx)
    gx = zeros(N,1);
end
if isempty(gy)
    gy = zeros(N,1);
end
if isempty(gz)
    gz = zeros(N,1);
end

% ----- calc b-value -----
if N>0
    g = [gx'; gy'; gz'];  % [Hz/m = 1/m/s]
    q = cumsum(g,2) * dt; % [Hz/m*s = 1/m]
    Q = 0; % diffusion encoding tensor
    for j=1:size(q,2)
        Q = Q + q(:,j) * q(:,j)' * dt; % [s/m^2]
    end
    Q = Q / 1e6; % [s/mm^2]
    b = trace(Q); % b-value of the diffusion encoding (isotropic case)
else
    b = 0;
end

end