%% init pulseq
clear
seq_name = 'mrf_';

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
FOV.Nxy    = 256;        % [ ] matrix size
FOV.Nz     = 1;          % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy = 256 *1e-3;  % [m] FOV geometry
FOV.dz     = 5   *1e-3;  % [m] slice thickness
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
SPI.mode_2D_3D     = '2D';      % '2D', '3D' or '3D_stacked'
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
SPI.spoil_nTwist   = 4;          % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode  = 'lin';      % 'lin' or 'quad'
SPI.spoil_rf_inc   = 0 *pi/180;  % [rad] rf phase increment
SPI.spoil_duration = 1.0 *1e-3;  % [s] slice spoiler duration

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

%% params: Inversion
INV.rf_type      = 'HYPSEC_inversion';
INV.tExc         = 10 *1e-3;  % [s]  hypsech pulse duration
INV.beta         = 700;       % [Hz] maximum rf peak amplitude
INV.mu           = 4.9;       % [ ]  determines amplitude of frequency sweep
INV.inv_rec_time = 0.01;      % [s]  inversion recovery time
INV = INV_init(INV, FOV, system);

%% noise pre-scans
SPI.Nnoise = 16;
SPI_add_prescans();

%% create sequence

% inversion
seq.addTRID('inversion');
INV_add();

% spiral imaging
for loop_NR = 1:SPI.NR
    SPI_add();
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();