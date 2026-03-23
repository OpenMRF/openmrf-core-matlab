%% load mrf study
clear

% enter paths of mrf and trajectory rawdata
% note: GE, United Imaging and Philips need explicit description of the
% pulseq backup paths. for Siemens, rawdata and pulseq backups can be
% automatically linked

vendor           = 'Siemens';
path_raw_mrf     = 'Q:/data/Pulseq/Rawdata/mgram/SolaMIITT/250311_tests_v151/meas_MID00308_FID84185_260311_0123_mgram_mrf_yun.dat';
path_raw_traj    = 'Q:/data/Pulseq/Rawdata/mgram/SolaMIITT/250311_tests_v151/meas_MID00321_FID84198_260311_0203_mgram_traj_260311_0123.dat';
path_backup_mrf  = [];
path_backup_traj = [];

%% load MRF rawdata
[DATA, NOISE, PULSEQ, study_info] = pulseq_read_meas(path_raw_mrf, path_backup_mrf, vendor);

%% load measured trajectory
if ~isempty(path_raw_traj)
    [ktraj_meas, ktraj_hash] = TRAJ_reco(path_raw_traj, path_backup_traj, vendor, [8:8:48]);
    if ~strcmp(ktraj_hash, pulseq_get_wave_hash(SPI_load_ktraj(PULSEQ)))
        ktraj_meas = []; % ensure compatibility of trajectory scan
    end
else
    ktraj_meas = [];
end

%% parameters: dictionary simulation
params_dict.seq_path    = [];       % full path of .seq file; [] for automatic search in user backups
params_dict.sim_mode    = 'BLOCH';  % 'BLOCH' or 'EPG'
params_dict.comp_energy = 0.9999;   % svd compression energy: 0 for uncompressed dictionary
params_dict.N_iso       = 1000;     % number of isochromats for bloch simulation
params_dict.s_fac       = 2;        % factor for out-of-slice simulation
params_dict.f0          = [];       % larmor frequency f0: required to simulate rf pulses with ppm freq offsets
params_dict.time_stamps = study_info.time_stamps; % adc times stamps: required for correction of trigger delays in cMRF
params_dict.soft_delays = [];       % soft delay user input: required for correction of the acq window in cMRF
params_dict.flag_kz     = [];       % find kz partitions for stacked 3D MRF -> eliminate unnecessary partitions
params_dict.echo_mode   = [];       % echo mode; default: 'spiral_out'

%% parameters: dictionary and look-up table

% T1, T2
P.T1.range = [0.01,  4]; P.T1.factor = 1.05;
P.T2.range = [0.001, 3]; P.T2.factor = 1.05;
P = MRF_get_param_dict(P, {'T2<T1'});
look_up       = [P.T1, P.T2];
look_up_names = {'T1', 'T2'};

% T1, T2, T1p
% P.T1.range  = [0.01,  4]; P.T1.factor  = 1.05;
% P.T2.range  = [0.001, 3]; P.T2.factor  = 1.05;
% P.T1p.range = [0.001, 3]; P.T1p.factor = 1.05;
% P = MRF_get_param_dict(P, {'T2<T1', 'T2<T1p', 'T1p<T1'});
% look_up       = [P.T1, P.T2, P.T1p];
% look_up_names = {'T1', 'T2', 'T1p'};

% T1, T2, B1+ correction
% P.T1.range  = [0.01,  4]; P.T1.factor = 1.025;
% P.T2.range  = [0.001, 3]; P.T2.factor = 1.025;
% P.db1.range = [0.8, 1.2]; P.db1.step  = 0.025;
% P = MRF_get_param_dict(P, {'T2<T1'});
% look_up       = [P.T1, P.T2, P.db1];
% look_up_names = {'T1', 'T2', 'db1'};

params_dict.P             = P;
params_dict.look_up       = look_up;
params_dict.look_up_names = look_up_names;
clear P look_up look_up_names;

%% parameters: low-rank reconstruction & matching
params_reco.TestReco           = true;   % do test reco: mixed contrast of all TRs
params_reco.DirectMatching     = false;  % do reco via direct matching
params_reco.DirectMatching_SVD = false;  % do reco via SVD compression of the dictionary before matching
params_reco.LowRank            = true;   % do reco via iterative low-rank reconstruction
params_reco.ROVIR              = true;   % use ROVIR coils for outer FOV artifact suppression
params_reco.rovir_mode         = 'auto'; % 'auto' uses the oversampled region; 'manual' interactive ROI selection
params_reco.rovir_method       = 'auto'; % 'auto' automatic estimation of virtual coils based on thresholding; 'manual' define number of coils
params_reco.rovir_thresh       = 2;      % threshold for automatic coil compression
params_reco.CoilComp           = true;   % additional SVD coil compression after ROVIR
params_reco.NCoils_v           = 8;      % final number of virtual coils after ROVIR and SVD compression
params_reco.ESPIRiT            = true;   % use ESPIRiT or openadapt for calculating cmaps
params_reco.readOS             = 2;      % read oversampling factor
params_reco.NBlocks            = 12;     % number of blocks for pattern matching
params_reco.zero_params.on_off = true;   % parameters for zero interpolation filling
params_reco.zero_params.filter = 'kaiser';
params_reco.zero_params.factor = 2.0;
params_reco.zero_params.radius = 6.0;

% low-rank reconstruction parameters
params_LR = setupParameters_LowrankMRF2D();

%% simulate dictionary, reconstruction of rawdata & match parameters
[match, images] = mrf_reco(DATA, NOISE, PULSEQ, params_dict, params_reco, params_LR, ktraj_meas);

%% visualize match results
vis = 'LR'; % 'direct' | 'SVD' | 'LR'
matchVis = match.(vis);

% --- limits + colormaps ---
t1lims  = [0 2000] * 1e-3;
t2lims  = [0 1000] * 1e-3;
t1plims = [0 1000] * 1e-3;
db1lims = 0.2 * [-1, 1] + 1;
t1cmp   = get_cmp('T1', 1000, 1);
t2cmp   = get_cmp('T2', 1000, 1);
t1pcmp  = get_cmp('T2', 1000, 1);
db1cmp  = get_cmp('blue_red', 1000);    

% Figure 1: M0 + T1 + T2
figure();
ax1 = subplot(1,3,1);
imagesc(abs(matchVis.M0)); axis image off; colormap(ax1, gray); colorbar;
title(sprintf('M0 %s', vis), 'Interpreter','none');

ax2 = subplot(1,3,2);
imagesc(matchVis.T1, t1lims); axis image off; colormap(ax2, t1cmp); colorbar;
title(sprintf('T1 %s', vis), 'Interpreter','none');

ax3 = subplot(1,3,3);
imagesc(matchVis.T2, t2lims); axis image off; colormap(ax3, t2cmp); colorbar;
title(sprintf('T2 %s', vis), 'Interpreter','none');

linkaxes([ax1 ax2 ax3]);
clear ax1 ax2 ax3;

% optional Figure 2: T1p (only if present)
if isfield(matchVis, 'T1p') && ~isempty(matchVis.T1p)
    figure();
    imagesc(matchVis.T1p, t1plims); axis image off; colormap(gca, t1pcmp); colorbar;
    title(sprintf('T1p %s', vis), 'Interpreter','none');
end

% optional Figure 3: dB1+ (only if present)
if isfield(matchVis, 'db1') && ~isempty(matchVis.db1)
    figure();
    imagesc(matchVis.db1, db1lims); axis image off; colormap(gca, db1cmp); colorbar;
    title(sprintf('dB1+ %s', vis), 'Interpreter','none');
end

%% save results
res.match = match; res.images = images;
params.params_dict = params_dict; params.params_reco = params_reco; params.params_LR = params_LR;
save_study_results(res, params, study_info);