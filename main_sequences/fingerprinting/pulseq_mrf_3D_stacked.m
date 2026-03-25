%% init pulseq
clear
seq_name = 'mrf_3D_stacked';

% main flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation

% optional: select scanner
% pulseq_scanner = 'Siemens_Sola_1,5T_MIITT';

% optional: select pns sim orientation
% pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy    = 256;        % [ ] matrix size
FOV.Nz     = 48;         % [ ] numer of "stack-of-spirals"
FOV.fov_xy = 256 *1e-3;  % [m] FOV geometry
FOV.fov_z  = 96  *1e-3;  % [m] FOV geometry
FOV.dz     = FOV.fov_z * 0.8;  % [m] slab thickness, reduce to prevent aliasing along z
FOV_init();

%% params: MRF flipangles and repetition times
MRF.pattern         = 'yun'; % select pattern: e.g. 'yun' or 'cao'
[MRF.FAs, MRF.TRs]  = MRF_get_FAs_TRs(MRF.pattern, 1);
MRF.FAs(MRF.FAs==0) = 1e-6;
seq_name            = [seq_name MRF.pattern];
MRF.NR              = numel(MRF.FAs);

% or define a custom pattern
% MRF.pattern = 'sin_70';
% MRF.FAs     = MRF_calc_FAs_sin([5, 30, 200; 1, 70, 200; 10, 10, 100]) *pi/180;
% MRF.NR      = numel(MRF.FAs);
% MRF.TRs     = 12 *1e-3 *ones(MRF.NR,1);
% seq_name    = [seq_name MRF.pattern];

%% params: Spiral Readouts

% basic/import params for MRF
SPI.mode_2D_3D     = '3D_stacked';      % '2D', '3D' or '3D_stacked'
SPI.adcBW          = 400 *1e3;  % [Hz] desired receiver bandwidth
SPI.NR             = MRF.NR;    % [ ] number of repetitions
SPI.mrf_import.TRs = MRF.TRs;   % [s] repetition times
SPI.mrf_import.FAs = MRF.FAs;   % [rad] flip angles

% slice excitation
SPI.exc_mode      = 'sinc';      % 'sinc' or 'sigpy_SLR'
SPI.exc_shape     = 'ex';        % only for sigpy: 'st' or 'ex' 
SPI.exc_time      = 2.0 *1e-3;   % [s] excitation time
SPI.exc_tbw       = 6;           % [ ] time bandwidth product
SPI.exc_fa_mode   = 'import';    % 'equal',  'ramped',  'import'  
SPI.reph_duration = 1.0 *1e-3;   % [s] slice rephaser duration

% params: gradient spoiling & rf spoiling
SPI.spoil_nTwist   = 2 * FOV.Nz; % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode  = 'lin';      % 'lin' or 'quad'
SPI.spoil_rf_inc   = 0 *pi/180;  % [rad] rf phase increment
SPI.spoil_duration = 1.2 *1e-3;  % [s] slice spoiler duration

% k-space geometry params
SPI.geo.design        = 'spiral';
SPI.geo.design_fun    = 'hargreaves';
SPI.geo.N_interleaves = 24;
SPI.geo.FOV_coeff     = [1 -0.5];
SPI.geo.Ns            = 1e5;
SPI.geo.ds            = 1e-3;
SPI.geo.flag_rv       = 0;

% k-space projection params
SPI.proj.mode = 'RoundGoldenAngle'; % 'Equal2Pi' 'RoundGoldenAngle'
SPI.proj.Nid  = 48; % number of unique projections

% gradient & slew rate limitation factors
SPI.lim_exc_slew   = 0.9; % slice excitation
SPI.lim_reph_grad  = 0.9; % slice rephasing
SPI.lim_reph_slew  = 0.6; % slice rephasing
SPI.lim_read_grad  = 0.9; % readout
SPI.lim_read_slew  = 0.5; % readout
SPI.lim_spoil_grad = 0.9; % slice spoiling
SPI.lim_spoil_slew = 0.6; % slice spoiling

[SPI, ktraj_ref] = SPI_init(SPI, FOV, system, 1);
MRF.TRs = SPI.TR;

%% averages of center partitions
SPI.nz_zero  = find(SPI.kz_area==0);
SPI.kz_table = [ 1:SPI.nz_zero-4, ...
                 SPI.nz_zero-3, ...
                 SPI.nz_zero-2, SPI.nz_zero-2, ...
                 SPI.nz_zero-1, SPI.nz_zero-1, SPI.nz_zero-1, ...
                 SPI.nz_zero,   SPI.nz_zero,   SPI.nz_zero,   SPI.nz_zero, ...
                 SPI.nz_zero+1, SPI.nz_zero+1, SPI.nz_zero+1, ...
                 SPI.nz_zero+2, SPI.nz_zero+2, ...
                 SPI.nz_zero+3, ...
                 SPI.nz_zero+4:FOV.Nz];

rng('default');
[~, SPI.kz_rand] = sort(rand(size(SPI.kz_table)));
SPI.kz_table   = SPI.kz_table(SPI.kz_rand);

%% params: Saturation (magnetization reset)
SAT.mode            = 'on';
SAT.rf_type         = 'adiabatic_BIR4';
SAT.bir4_tau        = 10 *1e-3;  % [s]  bir4 pulse duration
SAT.bir4_f1         = 640;       % [Hz] maximum rf peak amplitude
SAT.bir4_beta       = 10;        % [ ]  am waveform parameter
SAT.bir4_kappa      = atan(10);  % [ ]  fm waveform parameter
SAT.bir4_dw0        = 30000;     % [rad/s] fm waveform scaling
SAT.sat_rec_time    = 3.0;       % [s] saturation recovery delay
SAT.crush_nTwists_z = FOV.Nz * 3.7;
SAT = SAT_init(SAT, FOV, system);

%% params: Inversion
INV.rf_type         = 'HYPSEC_inversion';
INV.tExc            = 10 *1e-3;  % [s]  hypsech pulse duration
INV.beta            = 700;       % [Hz] maximum rf peak amplitude
INV.mu              = 4.9;       % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time    = [10] *1e-3;
INV.crush_nTwists_z = FOV.Nz * 2.9;
INV = INV_init(INV, FOV, system);

%% noise pre-scans
SPI.Nnoise = 64;
SPI_add_prescans();

%% create sequence

% init START STOP labels
% used for dictionary calculation
mr.addCustomLabel('START'); % we simulate after the START label
mr.addCustomLabel('STOP'); % and the simulation ends at the STOP label
start_stop = 1;

% pre-saturation
SAT_add();

for loop_kz = SPI.kz_table
% for loop_kz = SPI.nz_zero % use this for tests with only the center kz partition

    % start simulation of .seq file (kz=0)
    if (loop_kz == SPI.nz_zero) && (start_stop==1)
        seq.addBlock(mr.makeLabel('SET', 'START', 1));
    end

    % saturation recovery
    SAT_add();

    % inversion
    INV_add();
    
    % spiral imaging
    for loop_NR = 1:SPI.NR
        SPI_add();
    end
    
    % stop simulation of .seq file (kz=0)
    if (loop_kz == SPI.nz_zero) && (start_stop==1)
        seq.addBlock(mr.makeLabel('SET', 'STOP', 1));
        start_stop = 0;
    end

end

%% plot sequence diagram
seq.plot('timeRange', [5 26]);

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();