%% init pulseq
% basic TSE (turbo-spin-echo) readout
clear
seq_name = 'tse_centric';

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
FOV.Nx       = 256;
FOV.Ny       = 256;
FOV.fov_x    = 256 *1e-3;
FOV.fov_y    = 256 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% TSE sequence parameters
TSE.Nrep      = 1;
TSE.n_echo    = 8;
TSE.TE        = 12    *1e-3;
TSE.TR        = 1000  *1e-3;
TSE.Ndummy    = 4;
TSE.exc_time  = 3.0 *1e-3;
TSE.rfc_time  = 3.0 *1e-3;
TSE.exc_tbw   = 4;
TSE.rfc_tbw   = 4;
TSE.t_acq     = 6.4 *1e-3 + 2*system.adcDeadTime;
TSE.os_mode   = 1;  % read oversampling: 0 off, 1 on
TSE.mode_exc  = 'sinc'; % 'sigpy_SLR' or 'sinc'
TSE.mode_rfc  = 'sinc'; % 'sigpy_SLR' or 'sinc'
TSE.enc_mode  = 'centric';

[TSE, ktraj_adc, ktraj_full] = TSE_init(TSE, FOV, system, 1);

%% add sequence blocks
ndummy = TSE.Ndummy;

for loop_rep = 1 : TSE.Nrep
for loop_TR  = 1-ndummy : TSE.nex
    TSE_add();
    ndummy = 0;
end
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();