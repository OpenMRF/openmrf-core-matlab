%% init pulseq
clear
seq_name = 'cmrf_t1_t2';

% main flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation

% optional: select scanner
pulseq_scanner = 'Siemens_Sola_1,5T_MIITT';

% optional: select pns sim orientation
% pns_orientation = 'all';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy    = 192;         % [ ] matrix size
FOV.fov_xy = 320  *1e-3;  % [m] FOV geometry
FOV.dz     = 8   *1e-3;   % [m] slice thickness
FOV_init();

%% params: MRF contrast encoding

% encoding list with contrast preparations for different segments
% 'No_Prep'    ->  use for recovery of longitudinal magnetization
% 'Saturation' ->  use for T1 encoding
% 'Inversion'  ->  use for T1 encoding
% 'T2'         ->  use for T2 encoding
% 'SL'         ->  use for T1p encoding
% 'MLEV'       ->  use for T2 or T2p encoding

MRF.enc_list = {
'Inversion';
'No_Prep';
'T2';
'T2';
'Inversion';
'No_Prep';
'T2';
'T2';
'Inversion';
'No_Prep';
'T2';
'T2';
'Inversion';
'No_Prep';
'T2';
'T2'
};

MRF.n_segm = numel(MRF.enc_list);

%% params: MRF flipangles and repetition times
MRF.nr     = 48;                         % numer of readouts per hear beat
MRF.NR     = MRF.n_segm * MRF.nr;        % total number of readouts
MRF.TRs    = 0.0 *1e-3 *ones(MRF.NR,1);  % minimize TRs
MRF.FA_min = 4 *pi/180;                  % [rad] minimum flip angle
MRF.FA_max = 15 *pi/180;                 % [rad] minimum flip angle
MRF.FAs    = MRF_calc_FAs_sin_rand(MRF.FA_min, MRF.FA_max, MRF.nr, MRF.n_segm);

%% params: MRF segment timings and trigger options

MRF.mode_seg = 'trigger'; % 'trigger', 'constant', 'variable'

switch MRF.mode_seg
    case 'trigger' % -> cardiac MRF        
        % this will enable the cardiac trigger; segment timings will be synchronized; soft delays will be used for adjusting the acq window on the scanner

    case 'constant' % -> constant segment durations (use e.g. for phantom validations or abdominal scans)
        MRF.seg_duration = 1000 *1e-3; % [s] all segments will have the same total duration including readouts; this example would emulate a 1000ms RR interval

    case 'variable' % -> variable recovery times (use e.g. for abdominal scans)
        MRF.rec_times = rand(MRF.n_segm, 1); % [s] segments can have different durations; choose a specific recovery time for each segment
end

%% params: Spiral Readouts

% basic/import params for MRF
SPI.mode_2D_3D     = '2D';      % '2D', '3D' or '3D_stacked'
SPI.adcBW          = 400 *1e3;  % [Hz] desired receiver bandwidth
SPI.NR             = MRF.NR;    % [ ] number of repetitions
SPI.nr             = MRF.nr;    % [ ] numer of readouts per hear beat
SPI.mrf_import.TRs = MRF.TRs;   % [s] repetition times
SPI.mrf_import.FAs = MRF.FAs;   % [rad] flip angles

% slice excitation
SPI.exc_mode      = 'sinc';      % 'sinc' or 'sigpy_SLR'
SPI.exc_shape     = 'ex';        % only for sigpy: 'st' or 'ex' 
SPI.exc_time      = 0.5 *1e-3;   % [s] excitation time
SPI.exc_tbw       = 2;           % [ ] time bandwidth product
SPI.exc_fa_mode   = 'import';    % 'equal',  'ramped',  'import'  
SPI.reph_duration = 0.3 *1e-3;   % [s] slice rephaser duration

% params: gradient spoiling & rf spoiling
SPI.spoil_nTwist   = 4 * FOV.Nz;  % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode  = 'lin';       % 'lin' or 'quad'
SPI.spoil_rf_inc   = 0 *pi/180;   % [rad] rf phase increment
SPI.spoil_duration = 0.5 *1e-3;   % [s] slice spoiler duration

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
SPI.lim_exc_slew   = 0.95; % slice excitation
SPI.lim_reph_grad  = 0.95; % slice rephasing
SPI.lim_reph_slew  = 0.95; % slice rephasing
SPI.lim_read_grad  = 0.90; % readout
SPI.lim_read_slew  = 0.55; % readout
SPI.lim_spoil_grad = 0.95; % slice spoiling
SPI.lim_spoil_slew = 0.95; % slice spoiling

[SPI, ktraj_ref] = SPI_init(SPI, FOV, system, 1);
MRF.TRs = SPI.TR;

%% params: Inversion
INV.rf_type      = 'HYPSEC_inversion';
INV.tExc         = 10 *1e-3;  % [s]  hypsech pulse duration
INV.beta         = 700;       % [Hz] maximum rf peak amplitude
INV.mu           = 4.9;       % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time = [0.01 35 380 130] *1e-3;
INV = INV_init(INV, FOV, system);

%% params: T2 preparation
T2.exc_mode   = 'adiabatic_BIR4';
T2.rfc_dur    = 2 *1e-3;   % [s]  duration of composite refocusing pulses
T2.bir4_tau   = 10 *1e-3;  % [s]  bir4 pulse duration
T2.bir4_f1    = 640;       % [Hz] maximum rf peak amplitude
T2.bir4_beta  = 10;        % [ ]  am waveform parameter
T2.bir4_kappa = atan(10);  % [ ]  fm waveform parameter
T2.bir4_dw0   = 30000;     % [rad/s] fm waveform scaling
T2.prep_times = [40 80 40 80 40 80 40 80] * 1e-3;  % [s] inversion times
T2            = T2_init(T2, FOV, system);

%% params: Fat Saturation
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

%% check MRF encoding params and adjust segment delays
MRF_check_enc_list();
MRF_adjust_segment_delays();

%% noise pre-scans
SPI.Nnoise = 16;
SPI_add_prescans();

%% create sequence
MRF_add_segments();

%% plot sequence diagram
seq.plot();

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();