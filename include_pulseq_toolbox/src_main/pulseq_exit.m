% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% set sequence definitions
seq.setDefinition('Name',    seq_name);
seq.setDefinition('Scan_ID', wip_id);

if isfield(FOV,'fov_x') && isfield(FOV,'fov_y') && isfield(FOV,'fov_z')
    seq.setDefinition('FOV', [FOV.fov_x FOV.fov_y FOV.fov_z]);
elseif isfield(FOV,'fov_xy') && isfield(FOV,'fov_z')
    seq.setDefinition('FOV', [FOV.fov_xy FOV.fov_xy FOV.fov_z]);
end

seq.setDefinition('Rot_Matrix', FOV.Rot_Mat);

%% calculate total sequence duration
TotalDuration = sum(seq.blockDurations);
disp(' ')
if TotalDuration > 180
    disp(['   acq time:  ' num2str(round(TotalDuration/60,1)) 'min'])
else
    disp(['   acq time:  ' num2str(round(TotalDuration,1)) 's'])
end

%% calculate maximum RF peak
RFpeak = zeros(numel(seq.rfLibrary.data),1);
for j=1:numel(seq.rfLibrary.data)
    RFpeak(j) = seq.rfLibrary.data(j).array(1);
end
RFpeak = max(RFpeak);
disp('  ')
disp(['   Maximum RF peak magnitude: ' num2str(round(RFpeak)) 'Hz  ->  ' num2str(round(RFpeak*2*pi/(2.67522*1e8)*1e6,2)) 'uT']);
if RFpeak>700
    disp('   high RF peak detected! Sequence might not work in strict mode!')
end
disp('  ')

%% check sequence timings and report
if flag_report==0
    timings_ok   = [];
    error_report = [];
    test_report  = [];
    warning('   no timing check performed!!!')
end

% timing check
if flag_report>0
    [timings_ok, error_report] = seq.checkTiming;
    disp(' ');
    if (timings_ok)
        disp('   Timing check passed successfully');
    else
        disp(   'Timing check failed! Error listing follows:');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    test_report = [];
end

% try to calculate grad spectrum test
if flag_report>0
    % try
    %     [gradSpec.R, gradSpec.Rax, gradSpec.F] = seq.gradSpectrum(system.ascfile);
    % end
end
if ~exist('gradSpec', 'var')
    gradSpec = [];
end

% full test report
if flag_report>1
    disp(' ');
    test_report = seq.testReport;
    for j=1:numel(test_report)
        fprintf(test_report{1,j});
    end
    clear j
end

%% calculate PNS stimulation
if ~exist('flag_pns', 'var')
    flag_pns = 0;
end
if flag_pns==1
    if isempty(system.ascfile)
        warning('no .asc file for PNS simulation found!');
        flag_pns = 0;
    else
        if ~exist('pns_orientation', 'var')
            [pns_ok, pns_norm, pns_comp] = PNS_sim(seq, 'axial');
        else
            [pns_ok, pns_norm, pns_comp] = PNS_sim(seq, pns_orientation);
        end
    end
else
    pns_ok   = [];
    pns_norm = [];
    pns_comp = [];
end

%% calculate gradient acoustic spectrum
if ~exist('flag_sound', 'var')
    flag_sound = 0;
end
if flag_sound==1
    seq.sound();
end

%% write .seq file and save some backup files
if flag_backup>0

    % display sequence name for export
    disp(' ');
    disp(['   seq name for export:  ' seq_name]);
    
    % create paths for export and backup
    backup_path = [ pulseq_path '/Pulseq_Workspace/' pulseq_user '/' seq_id(1:6) ];
    if ~exist(backup_path, 'dir') && flag_backup>0
       mkdir(backup_path);
    end
    backup_path   = [backup_path '/' seq_id];
    external_path = [ backup_path '/' seq_name '.seq' ];

    % check for duplicate pulseq backup
    check_duplicate = isfile([ backup_path '/backup_' seq_id '_workspace.mat' ]);    
    if check_duplicate
        warning('pulseq duplicate detected! .seq writing stopped!');
    else

        % write .seq file
        mkdir(backup_path);
        seq.write(external_path);
        
        % backup: matlab code
        try
            copyfile( filepath, [ backup_path '/backup_' seq_name '.m' ] );
        end
        % backup: header information
        PULSEQ.seq_name       = seq_name;
        PULSEQ.seq_id         = seq_id;
        PULSEQ.scan_id        = scan_id;
        PULSEQ.md5_hash       = seq.signatureValue;
        PULSEQ.pulseq_lab     = pulseq_lab;
        PULSEQ.pulseq_scanner = pulseq_scanner;
        PULSEQ.pulseq_user    = pulseq_user;
        PULSEQ.pulseq_path    = pulseq_path;
        PULSEQ.backup_path    = backup_path;
        PULSEQ.seq            = seq;
        PULSEQ.system         = system;
        PULSEQ.gradSpec       = gradSpec;
        PULSEQ.timings_ok     = timings_ok;
        PULSEQ.error_report   = error_report;
        PULSEQ.test_report    = test_report;
        PULSEQ.duration       = TotalDuration;
        PULSEQ.FOV            = FOV;
        PULSEQ.RFpeak         = RFpeak;

        % backup: scan specific objects and variables
        backup_list = {
        'ADIASL';
        'AHP';
		'BSSFP';
        'EPI';
        'FAT';
        'GRE';
        'INV';
        'MLEV';
        'MRF';
        'PRESS';
        'RAD';
        'REX';
        'SAT';
        'SL';
        'SPI';
        'SPITSE';
        'TSE';
        'TRAJ';
        'TREX';
        'TRIG_IN';
        'TRIG_OUT';
        'T2';
        'UTE';
        'WASABI';
        'ktraj_ref';
        'ktraj_reco';
        'ktraj_adc';
        'ktraj_full'
        };
        for j=1:numel(backup_list)
            eval(['if exist(''' backup_list{j} ''', ''var''); PULSEQ.' backup_list{j} ' = ' backup_list{j} '; end;']);
        end
        save( [ backup_path '/backup_' seq_id '_workspace.mat' ], 'PULSEQ', '-v7.3');
        
        % backup toolbox functions
        git_path = pulseq_get_path('pulseq_init');
        git_path = git_path(1:end-9);
        tar([backup_path '/include_pulseq_toolbox.tar'], git_path);

        % write txt file with md5 identification and header
        writelines( {PULSEQ.md5_hash; ...
                     PULSEQ.pulseq_lab; ...
                     PULSEQ.pulseq_scanner; ...
                     PULSEQ.seq_name; ...                     
                     PULSEQ.backup_path; ...
                     [num2str(PULSEQ.duration,'%.3f') 's,  ' num2str(PULSEQ.duration/60,'%.1f') 'min']}, ...
                    [backup_path '/md5_' PULSEQ.md5_hash '.txt']);
        
        % copy .seq file to Seq_Export folder
        if (flag_backup>1)
            if ~exist([pulseq_path '/Seq_Export/' pulseq_user], 'dir')
               mkdir([pulseq_path '/Seq_Export/' pulseq_user]);
            end
            copyfile( external_path, [pulseq_path '/Seq_Export/' pulseq_user '/' seq_name '.seq'] );
        end
        
        % send .seq file to hydra -> only for users from EP5 lab
        if (flag_backup>2)
            pulseq_scp_hydra(pulseq_user, external_path);
        end

    end
end

%% export RF adjustment scan for GE
if flag_GE == 1 && flag_backup>0
    if exist('SPI', 'var')
        GE_adj_receive_gain(system, 5, 2, SPI.adc, pi/2, FOV.dz, external_path, wip_id);
    elseif exist('GRE', 'var')
        GE_adj_receive_gain(system, 5, 2, GRE.adc, pi/2, FOV.dz, external_path, wip_id);
    elseif exist('TRAJ', 'var')
        GE_adj_receive_gain(system, 5, 2, TRAJ.adc, pi/2, 0.01, external_path, wip_id);
    end
end

%% fast EPG simulation of sequence
if flag_mrf == 1
    if ~exist('external_path', 'var')
        external_path = [ pulseq_path '/Pulseq_Workspace/' pulseq_user '/' seq_id(1:6) ];
        external_path = [external_path '/' seq_id];
        external_path = [ external_path '/' seq_name '.seq' ];
    end
    MRF_sim_on_the_fly(seq, external_path, flag_mrf);
end
