function [rawdata, noise, PULSEQ, study_info] = pulseq_read_meas(path_raw, path_backup, vendor)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 22.03.2026
    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V3, 25.08.2026
    
    % ----- Input: -----
    % path_raw:      file path to rawdata (.dat for Siemens, .h5 for ISMRMRD, .mat for others)
    % path_backup:   file path to PULSEQ backup (.mat, not necessary for Siemens)
    % vendor:        'Siemens'  'UI'  'GE'  'Philips'  'ISMRMRD'
    
    % ----- Output: -----
    % rawdata:     [coils x tr x adc] complex rawdata
    % noise:       [coils x rep x adc] complex noise pre-scans
    % PULSEQ       struct containing all sequence parameters
    % study_info:  struct containing scan/study specific header information

    %% init paths

    % defaults
    if nargin<3
        vendor = [];
    end
    if nargin<2
        path_backup = [];
    end
    if nargin<1
        path_raw = [];
    end

    % select rawdata via uigetfile()
    if isempty(path_raw)
        [~, ~, path_raw]       = pulseq_get_user_definitions([], 0);
        [temp_name, temp_path] = uigetfile('*.*', 'Select a file', path_raw );
        path_raw               = [temp_path temp_name];
        clear temp_name temp_path;
    end

    % check file extensions & vendors
    [~,~,temp_ext] = fileparts(path_raw);
    if isempty(temp_ext)
        error('path_raw needs a file extension .dat .h5 or .mat');
    else
        switch temp_ext
            case '.dat'
                vendor = 'Siemens';
            case '.h5'
                vendor = 'ISMRMRD';
            case '.mat'
                if isempty(vendor)
                    error('specify vendor for .mat rawdata');
                end
            otherwise
                error('rawdata has to be .dat .h5 or .mat');
        end
    end

    %% load rawdata and pulseq backups depending on vendor
    switch vendor
        case 'Siemens' % Siemens Healthcare
            [twix_obj, study_info, PULSEQ] = pulseq_read_meas_Siemens(path_raw);
            rawdata = permute(twix_obj.image.unsorted(), [2, 3, 1]); % coils x tr x adc
            if ~isempty(path_backup)
                load(path_backup);
                warning('automatic PULSEQ backup was overwritten for Siemens scan!');
            end            

        case 'UI' % United Imaging Healthcare
            [rawdata, study_info, PULSEQ] = pulseq_read_meas_UI(path_raw);
            rawdata = (rawdata(:,:,1:2:end) + rawdata(:,:,2:2:end)) / 2; % remove oversampling
            if ~isempty(path_backup)
                load(path_backup);
                warning('automatic PULSEQ backup was overwritten for United Imaging scan!');
            end

        case 'GE' % General Electric
            [rawdata, study_info, PULSEQ] = pulseq_read_meas_GE(path_raw);
            if ~isempty(path_backup)
                load(path_backup);
                warning('automatic PULSEQ backup was overwritten for GE scan!');
            end

        case 'Philips' % Koninklijke Philips
            [rawdata, study_info, PULSEQ] = pulseq_read_meas_Philips(path_raw); % to do
            if ~isempty(path_backup)
                load(path_backup);
                warning('automatic PULSEQ backup was overwritten for Philips scan!');
            end    

        case 'ISMRMRD' % ISMRMRD
            [rawdata, study_info, PULSEQ] = pulseq_read_meas_ISMRMRD(path_raw);
            if ~isempty(path_backup)
                load(path_backup);
                warning('automatic PULSEQ backup was overwritten for ISMRMRD scan!');
            end
    end

    %% split meas data and noise pre-scans
    if isfield(PULSEQ, 'SPI') && isfield(PULSEQ.SPI, 'Nnoise') && PULSEQ.SPI.Nnoise>0
        noise = rawdata(:,1:PULSEQ.SPI.Nnoise,:);
        rawdata(:,1:PULSEQ.SPI.Nnoise,:) = [];
    else
        noise = [];
    end

    %% check adc time stamps and soft delays
    if ~isfield(study_info, 'time_stamps')
        study_info.time_stamps = [];
    end
    if ~isfield(study_info, 'soft_delays')
        study_info.soft_delays = [];
    end
    if ~isempty(study_info.time_stamps) % if noise prescans have no time stamps
        if numel(study_info.time_stamps) < (size(noise,2)+size(rawdata,2))
            n_dummy = (size(noise,2)+size(rawdata,2)) - numel(study_info.time_stamps);
            study_info.time_stamps = [ (1:n_dummy)'*1e-12; study_info.time_stamps ]; % we fix it with nearly zero values
        end
    end

end
