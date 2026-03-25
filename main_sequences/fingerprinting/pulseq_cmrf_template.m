%% NOTE
% This script is a template for building a cardiac MRF sequence.
% It demonstrates how to combine different preparation modules and readouts.
% The specific choice of preparations here is not intended for meaningful experiments.
% For realistic implementations, please refer to:
%   - pulseq_cmrf_t1_t2
%   - pulseq_cmrf_t1_t2_t1p

%% init pulseq
clear
seq_name = 'cardiac_mrf';

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

MRF.enc_list = { % the following encoding list is only an example which includes all possible preparations
'Inversion';
'No_Prep';
'T2';
'T2';
'Inversion';
'Saturation';
'SL';
'SL';
'Inversion';
'Saturation';
'No_Prep';
'No_Prep';
'Inversion';
'No_Prep';
'MLEV';
'MLEV'
};

MRF.n_segm = numel(MRF.enc_list);

%% params: MRF flipangles and repetition times
MRF.nr     = 48;                         % numer of readouts per hear beat
MRF.NR     = MRF.n_segm * MRF.nr;        % total number of readouts
MRF.TRs    = 0.0 *1e-3 *ones(MRF.NR,1);  % minimize TRs
MRF.FA_min = 4 *pi/180;                  % [rad] minimum flip angle
MRF.FA_max = 15 *pi/180;                 % [rad] minimum flip angle
MRF.FAs    = MRF_calc_FAs_sin_rand(MRF.FA_min, MRF.FA_max, MRF.nr, MRF.n_segm);

% or use constant FAs
% MRF.FAs = 15*pi/180 * ones(MRF.NR,1);

% or import from .mat file
% load('FAs.mat');
% MRF.FAs = FAs;

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
SPI.spoil_nTwist   = 4;           % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode  = 'lin';       % 'lin' or 'quad'
SPI.spoil_rf_inc   = 0 *pi/180;   % [rad] rf phase increment
SPI.spoil_duration = 0.5 *1e-3;  % [s] slice spoiler duration

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
INV.inv_rec_time = [15 75 150 250] *1e-3;
INV = INV_init(INV, FOV, system);

%% params: Saturation
SAT.mode         = 'on';
SAT.rf_type      = 'adiabatic_BIR4';
SAT.bir4_tau     = 10 *1e-3;  % [s]  bir4 pulse duration
SAT.bir4_f1      = 640;       % [Hz] maximum rf peak amplitude
SAT.bir4_beta    = 10;        % [ ]  am waveform parameter
SAT.bir4_kappa   = atan(10);  % [ ]  fm waveform parameter
SAT.bir4_dw0     = 30000;     % [rad/s] fm waveform scaling
SAT.sat_rec_time = [100 200] *1e-3; % [s] saturation times
SAT = SAT_init(SAT, FOV, system);

%% params: T2 preparation
T2.exc_mode   = 'adiabatic_BIR4';
T2.rfc_dur    = 2 *1e-3;   % [s]  duration of composite refocusing pulses
T2.bir4_tau   = 10 *1e-3;  % [s]  bir4 pulse duration
T2.bir4_f1    = 640;       % [Hz] maximum rf peak amplitude
T2.bir4_beta  = 10;        % [ ]  am waveform parameter
T2.bir4_kappa = atan(10);  % [ ]  fm waveform parameter
T2.bir4_dw0   = 30000;     % [rad/s] fm waveform scaling
T2.prep_times = [40 80] * 1e-3;  % [s] inversion times
T2            = T2_init(T2, FOV, system);

%% params: Spin-Lock

% spin-lock pulses
SL.relax_type = {'T1p'};         % T1p or T2p or T2
SL.seq_type   = {'BSL'};         % BSL or CSL or RESL
SL.tSL        = [40 80] *1e-3;   % [s]  SL time
SL.fSL        = [200 200];       % [Hz] SL amplitude

% excitation pulses
SL.exc_mode  = 'adiabatic_AHP';  % 'adiabatic_AHP', 'sinc', 'sigpy_SLR' or 'bp'
SL.exc_time  = 3.0 *1e-3;        % [s] excitation time
SL.adia_wmax = 600 * 2*pi;       % [rad/s] amplitude of adiabatic pulse

% refocusing pulses
SL.rfc_mode = 'bp';              % 'bp', 'sinc', 'sigpy_SLR' or 'comp'
SL.rfc_time = 1.0 *1e-3;         % [s] refocusing time

SL = SL_init(SL, FOV, system);

%% params: MLEV T2p preparation
MLEV.n_mlev   = [4 8];           % number of MLEV4 preps
MLEV.fSL      = 250;             % [Hz] eff spin-lock field strength
MLEV.t_inter  = 1 *1e-5;         % [s]  inter pulse delay for T2 preparation
MLEV.exc_mode = 'adiabatic_AHP'; % 'adiabatic_BIR4' or 'adiabatic_AHP'
MLEV = MLEV_init(MLEV, FOV, system);

%% params: Fat Saturation
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

%% check MRF encoding params
MRF_check_enc_list();

%% adjust dynamic segment delays
MRF_adjust_segment_delays();

%% Trigger Mode
% on:  Cardiac
% off: Abdominal or Phantom

MRF.mode_trig = 'off';

% calc fixed segment timings
MRF.acq_duration      = sum(SPI.TR(1:MRF.nr));
MRF.prep_acq_duration = MRF.prep_max + MRF.acq_duration;

if strcmp(MRF.mode_trig, 'on')
    MRF.delay_soft = mr.makeSoftDelay(0, 'acq_end', 'offset', -MRF.prep_acq_duration, 'factor', 1); % block_duration [s] = offset [s] + input [s] / factor
    TRIG_IN = mr.makeTrigger('physio1', 'system', system, 'delay', 10e-6, 'duration', 10e-3); % Input Trigger
else
    MRF.seg_duration = 1000 *1e-3; % [s] adjust segment duration
    MRF.delay_soft   = mr.makeDelay( round((MRF.seg_duration-MRF.prep_acq_duration)/system.gradRasterTime)*system.gradRasterTime ); % fixed delay
end

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