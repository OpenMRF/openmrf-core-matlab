function [match, images, dict_comp] = mrf_reco(DATA, NOISE, PULSEQ, params_dict, params_reco, params_LR, ktraj_meas)
 
% ----- Input -----
% DATA:          [NCoils x NR x NRead]
% NOISE:         [NCoils x Nnoise x NRead]
% PULSEQ:        structured object with pulseq backup data
% params_dict:   structured object with dictionary calculation settings
% params_reco:   structured object with general reconstruction settings
% params_LR:     structured object for Low Rank settings
% ktraj_meas:    [2 x Nid x NRead] -> 1/m

% ----- Output -----
% match:     structured object with results of matching
% images:    structured object with reconstructed images
% dict_comp: structured object with compressed dictionary

% show mrf IP notice
mrf_ip_notice();

% only for cMRF: control motion with sliding window reconstruction
if isfield(PULSEQ.MRF, 'n_segm') && params_reco.DirectMatching
    params_reco.NSeg = PULSEQ.MRF.n_segm;
end

% default for zero interpolation filling
if ~isfield(params_reco, 'zero_params')
    params_reco.zero_params.on_off = false;
end

% save uncompressed image series from LR reconstruction
if ~isfield(params_LR, 'save_uncomp')
    params_LR.save_uncomp = false;
end  

%% read .seq file and generate simulation instructions

if isempty(params_dict.seq_path)
    [~, ~, params_dict.seq_path] = pulseq_get_user_definitions();
    params_dict.seq_path = [params_dict.seq_path '/Pulseq_Workspace/' PULSEQ.pulseq_user '/' PULSEQ.seq_id(1:6) '/' PULSEQ.seq_id '/' PULSEQ.seq_name '.seq' ];
end

[~, SIM] = MRF_read_seq_file( params_dict.seq_path, ...     % path of the original .seq file which was measured
                              params_dict.f0, ...           % larmor frequency f0; required to simulate rf pulses with ppm freq offsets
                              params_dict.time_stamps, ...  % adc times stamps; used for correction of trigger delays
                              params_dict.soft_delays,...   % soft delay input; used for correction of sequence timings
                              params_dict.flag_kz, ...      % find kz partitions for stacked 3D MRF -> eliminate unnecessary partitions
                              params_dict.echo_mode, ...    % echo mode; default: 'spiral_out'
                              PULSEQ.SPI.adcNPad,...        % adc padding
                              PULSEQ.SPI.adc.dwell,...      % adc dwell time
                              1e-6, ...                     % raster time for the simulation 
                              0);                           % flag_plot

%% caclulate dictionary
switch params_dict.sim_mode
    case 'EPG'
        SIM  = MRF_sim_pre(SIM, params_dict.P, [], 'EPG', 0, 0);
        dict = MRF_sim_EPG(SIM, params_dict.P, params_dict.comp_energy);
    case 'BLOCH'
        if ~isfield(params_dict, 'z')
            if ~isfield(params_dict, 's_fac')
                params_dict.s_fac = 2.0;
            end
            params_dict.z = linspace(-1/2, 1/2, params_dict.N_iso)' * params_dict.s_fac * PULSEQ.FOV.dz;
        end
        SIM  = MRF_sim_pre(SIM, params_dict.P, params_dict.z, 'BLOCH', 1, 0); 
        dict = MRF_sim_BLOCH(SIM, params_dict.P, params_dict.z, [], params_dict.comp_energy);
end

%% load k-space trajectory
fov = PULSEQ.FOV.fov_xy; % [1/m]
Nxy = PULSEQ.FOV.Nxy;   
if isempty(ktraj_meas)
    ktraj = SPI_load_ktraj(PULSEQ);
    if size(ktraj,1) == 3 && ~strcmp(PULSEQ.SPI.mode_2D_3D, '3D')
        ktraj(3,:,:) = [];
    end
else
    ktraj = ktraj_meas;
end
if isfield(PULSEQ.SPI, 'proj')
    proj_id = PULSEQ.SPI.proj.id;
else
    proj_id = PULSEQ.SPI.phi_id;
end
if isfield(params_reco, 'skope_full_traj') && params_reco.skope_full_traj
    proj_id = (1:size(ktraj,2))';
end
clear ktraj_meas;

%% noise pre-whitening
if ~isempty(NOISE)
    DATA = mg_noise_prewhitening(DATA, NOISE(:,:,20:end), 'cholesky', 1);
end
clear NOISE;

%% adc padding
if isfield(PULSEQ.SPI, 'adcNPad')
    if numel(PULSEQ.SPI.adcNPad) == 2
        DATA  = DATA(:,:, PULSEQ.SPI.adcNPad(1):PULSEQ.SPI.adcNPad(2));
        ktraj = ktraj(:,:, PULSEQ.SPI.adcNPad(1):PULSEQ.SPI.adcNPad(2));
    else
        DATA  = DATA(:,:,PULSEQ.SPI.adcNPad+1:end);
        ktraj = ktraj(:,:,PULSEQ.SPI.adcNPad+1:end);
    end
end

%% optional: try mixed test reco & visulaize
if params_reco.TestReco
    test_reco = mrf_mixed_reco(DATA, ktraj(:,proj_id,:), Nxy, fov);
end

%% reshape data; prepare reconstruction
DATA    = single(permute(DATA, [3,1,2]));
kx      = single(squeeze(ktraj(1,:,:)))';
ky      = single(squeeze(ktraj(2,:,:)))';
look_up = single(params_dict.look_up);
if ~isstruct(dict)
    dict = single(dict);
end
clear ktraj;

% init
match               = struct();
images              = struct();
dict_comp           = struct();
[NRead, NCoils, NR] = size(DATA);
Nid                 = size(kx,2);

% crop the center 1x FOV for display
center_pix = Nxy*params_reco.readOS/2 + (-Nxy/2 + 1 : Nxy/2);
params_LR.pixelRangeToShow = center_pix;

%% SVD Dictionary Compression
if ~isstruct(dict)
    % See IEEE TMI paper by Debra McGivney
    [~, S ,V] = svd(dict.', 'econ');
    if ~isfield(params_reco, 'NPCs')
        svals = diag(S)/S(1);
        svals = cumsum(svals.^2./sum(svals.^2));
        NPCs  = find(svals > 0.9999, 1, 'first');
    else
        NPCs  = params_reco.NPCs;
    end
    dict_phi       = single(V(:,1:NPCs));
    dict_svd       = (dict.'*dict_phi).';
    dict_comp.phi  = dict_phi;
    dict_comp.svd  = dict_svd;
    dict_comp.NPCs = NPCs; 
    clear S V svals;    
else
    dict_phi  = dict.phi;
    dict_svd  = dict.svd;
    NPCs      = dict.NPCs;
    dict_comp = dict;
    if params_reco.DirectMatching
        params_reco.DirectMatching = false;
        warning('direct matching disabled! not possible with compressed dictionary');
    end
end

%% Set up the NUFFT for gridding

% calculate global density compensation function
kmax       = max(sqrt(kx(:).^2 + ky(:).^2));
kx         = kx/kmax*Nxy/2;    
ky         = ky/kmax*Nxy/2;
Ni         = size(true(Nxy,Nxy));
nufft_args = {Ni, [6 6], 2*Ni, Ni/2, 'table', 2^12, 'minmax:kb'};
G          = Gmri([kx(:) ky(:)]/fov, true(Nxy,Nxy), 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
dcf_all    = abs(mri_density_comp([kx(:) ky(:)]/fov, 'pipe', 'fix_edge', 2, 'G', G.arg.Gnufft));
dcf_all    = reshape(dcf_all,size(kx));
clear kmax fov Ni nufft_args G

% calculate density compensation for each TR
ktraj = zeros(NRead, 1, NR, 'single');
dcf   = zeros(NRead, 1, NR, 'single');
for j = 1:NR
    ktraj(:,1,j) = kx(:, proj_id(j)) + 1i*ky(:, proj_id(j));
    dcf(:,1,j)   = dcf_all(:, proj_id(j));
end
ktraj = ktraj/Nxy;        % scale trajectory: -0.5 ... +0.5
dcf   = dcf/max(dcf(:));  % scale DCF: 0 ... 1

% built NUFFT operator
FT = NUFFT( kx/Nxy + 1i*ky/Nxy, ...
            dcf_all / max(dcf_all(:)), ...
            [0 0], ...
            [Nxy*params_reco.readOS Nxy*params_reco.readOS]);

%% ROVIR Outer FOV Artifact Suppression
if params_reco.ROVIR
    % See paper by Dauen Kim, "Region-optimized virtual (ROVir) coils", MRM 2021.

    % reco of 1st PC image before ROVir
    ksvd       = complex(zeros(NRead * Nid, NR, 'single'));
    images_pc1 = complex(zeros(Nxy * params_reco.readOS, Nxy * params_reco.readOS, NCoils, 'single'));
    for j = 1:NCoils
        for k=1:NR
            ksvd(NRead*(proj_id(k)-1) + 1:NRead*proj_id(k),k) = DATA(:,j,k);
        end
        kcomp      = reshape(ksvd*dict_phi, [NRead Nid NPCs]);
        image_temp = FT' * (kcomp(:,:,1) .* sqrt(dcf_all));
        images_pc1(:,:,j) = image_temp(:,:,1);
    end
    temp_before = abs(openadapt(permute(images_pc1, [3,1,2])));
    clear ksvd image_temp kcomp;

    if ~isfield(params_reco, 'rovir_mode')
        params_reco.rovir_mode = 'auto';
    end
    if strcmp(params_reco.rovir_mode, 'auto')
        % Here we use ROVIR to suppress signals from the oversampled region of the FOV (outside the 1x FOV)
        % Define ROI for signal within 1x FOV
        roiSignal = zeros(Nxy*params_reco.readOS, Nxy*params_reco.readOS);
        roiSignal(Nxy*params_reco.readOS/2 - Nxy/2+1 : Nxy*params_reco.readOS/2 + Nxy/2, ...
                  Nxy*params_reco.readOS/2 - Nxy/2+1 : Nxy*params_reco.readOS/2 + Nxy/2) = 1;
        % Define ROI for signal within outer 75% (25% transition band)
        buffer    = ceil(round(Nxy*1.25)/2)*2;
        roiInterf = zeros(Nxy*params_reco.readOS, Nxy*params_reco.readOS);
        roiInterf(Nxy*params_reco.readOS/2 - buffer/2 : Nxy*params_reco.readOS/2 + buffer/2-1, ...
                  Nxy*params_reco.readOS/2 - buffer/2 : Nxy*params_reco.readOS/2 + buffer/2-1) = 1;
        roiInterf = 1 - roiInterf;
        clear buffer;
    end
    if strcmp(params_reco.rovir_mode, 'manual')
        if ~isfield(params_reco, 'rovir_roiSignal') || ~isfield(params_reco, 'rovir_roiInterf')
            figure(679);
            hIm = imagesc(temp_before); axis image; axis off; colormap gray;
            title('Draw/adjust ellipse, then click Accept ROI');
            h   = drawellipse(gca);
            btn = uicontrol('Style', 'pushbutton', 'String', 'Accept ROI', ...
                            'Position', [20 20 100 30], ...
                            'Callback', @(src,evt) uiresume(gcbf));
                            uiwait(gcf);
                            drawnow;
            roiSignal = createMask(h, hIm) + 0;                
            buffer    = round(0.1 * size(images_pc1,1));
            roiInterf = (bwdist(roiSignal) > buffer) + 0;
            clear hIm h btn buffer;
            close(679);
        else
            roiSignal = params_reco.rovir_roiSignal;
            roiInterf = params_reco.rovir_roiInterf;
        end
    end

    % start ROVir coil beamforming
    switch params_reco.rovir_method
        case 'auto'
            DATA = ROVIR(DATA, images_pc1, params_reco.rovir_method, [], params_reco.rovir_thresh, true, roiSignal, roiInterf);
        case 'manual'
            DATA = ROVIR(DATA, images_pc1, params_reco.rovir_method, params_reco.rovir_NCoils_v, [], true, roiSignal, roiInterf);
    end
    NCoils = size(DATA, 2);
    
    % reco of 1st PC image after ROVir
    ksvd       = complex(zeros(NRead * Nid, NR, 'single'));
    images_pc1 = complex(zeros(Nxy * params_reco.readOS, Nxy * params_reco.readOS, NCoils, 'single'));
    for j = 1:NCoils
        for k=1:NR
            ksvd(NRead*(proj_id(k)-1) + 1:NRead*proj_id(k),k) = DATA(:,j,k);
        end
        kcomp      = reshape(ksvd*dict_phi, [NRead Nid NPCs]);
        image_temp = FT' * (kcomp(:,:,1) .* sqrt(dcf_all));
        images_pc1(:,:,j) = image_temp(:,:,1);
    end
    temp_after = abs(openadapt(permute(images_pc1, [3,1,2])));
    clear ksvd image_temp kcomp;

    % show ROIs before/after ROVir
    figure()
    subplot(2,3,1)
    imagesc(temp_before, [min(temp_before(:)) max(temp_before(:))]);
    axis image; axis off; colormap(gray(1000)); colorbar;
    title('Full Signal before ROVir')
    hold on
    contour(roiSignal, [0.5 0.5], 'g--', 'LineWidth', 1)
    contour(roiInterf, [0.5 0.5], 'r--', 'LineWidth', 1)
    subplot(2,3,2)
    imagesc(roiSignal .* temp_before, [0 max(temp_before(roiSignal==1))]);
    axis image; axis off; colormap(gray(1000)); colorbar;
    title('ROI Signal before ROVir')
    hold on
    contour(roiSignal, [0.5 0.5], 'g--', 'LineWidth', 1)
    subplot(2,3,3)
    imagesc(roiInterf .* temp_before, [0 max(temp_before(roiInterf==1))]);
    axis image; axis off; colormap(gray(1000)); colorbar;
    title('Interference Signal before ROVir')
    hold on
    contour(roiInterf, [0.5 0.5], 'r--', 'LineWidth', 1)
    subplot(2,3,4)
    imagesc(temp_after, [min(temp_after(:)) max(temp_after(:))]);
    axis image; axis off; colormap(gray(1000)); colorbar;
    title('Full Signal after ROVir')
    hold on
    contour(roiSignal, [0.5 0.5], 'g--', 'LineWidth', 1)
    contour(roiInterf, [0.5 0.5], 'r--', 'LineWidth', 1)
    subplot(2,3,5)
    imagesc(roiSignal .* temp_after, [0 max(temp_after(roiSignal==1))]);
    axis image; axis off; colormap(gray(1000)); colorbar;
    title('ROI Signal after ROVir')
    hold on
    contour(roiSignal, [0.5 0.5], 'g--', 'LineWidth', 1)
    subplot(2,3,6)
    imagesc(roiInterf .* temp_after, [0 max(temp_after(roiInterf==1))]);
    axis image; axis off; colormap(gray(1000)); colorbar;
    title('Interference Signal after ROVir')
    hold on
    contour(roiInterf, [0.5 0.5], 'r--', 'LineWidth', 1)
    clear roiSignal roiInterf images_pc1 temp_before temp_after;

end

%% SVD Coil Compression
if params_reco.CoilComp
    no_coils = min([params_reco.NCoils_v, NCoils]);
    [~, ~, v_coils] = svd(reshape(permute(DATA, [2,1,3]), NCoils, []).', 'econ');
    DATA     = reshape(permute(DATA, [1,3,2]), [NRead*NR NCoils]) * v_coils(:,1:no_coils);
    DATA     = permute(reshape(DATA, [NRead, NR, no_coils]), [1 3 2]);
    DATA     = single(DATA);
    disp(['SVD: keep ' num2str(no_coils) ' of ' num2str(NCoils) ' coils' ]);
    NCoils   = no_coils;
    clear no_coils v_coils;
end

%% Estimate coil sensitivity maps

% reco of 1st PC image
ksvd       = complex(zeros(NRead * Nid, NR, 'single'));
images_pc1 = complex(zeros(Nxy * params_reco.readOS, Nxy * params_reco.readOS, NCoils, 'single'));
for j = 1:NCoils
    for k=1:NR
        ksvd(NRead*(proj_id(k)-1) + 1:NRead*proj_id(k),k) = DATA(:,j,k);
    end
    kcomp      = reshape(ksvd*dict_phi, [NRead Nid NPCs]);
    image_temp = FT' * (kcomp(:,:,1) .* sqrt(dcf_all));
    images_pc1(:,:,j) = image_temp(:,:,1);
end
images.pc1 = images_pc1;
clear ksvd image_temp kcomp images_pc1;

% Estimate coil sensitivity maps from the 1st singular value MRF image
if params_reco.ESPIRiT
    % espirit coil sensitivities
    nyacs       = 24*params_reco.readOS;
    kcalib      = fft2c(images.pc1);
    kcalib      = kcalib(end/2-nyacs/2:end/2+nyacs/2-1,end/2-nyacs/2:end/2+nyacs/2-1,:);
    [kernel, s] = dat2Kernel(kcalib, [6 6]);
    eigThresh_k = 0.02;
    idx         = find(s >= s(1)*eigThresh_k, 1, 'last' );
    [M, ~]      = kernelEig(kernel(:,:,:,1:idx), [Nxy*params_reco.readOS Nxy*params_reco.readOS]);
    cmaps       = M(:,:,:,end);
    clear nyacs kcalib kernel s eigThresh_k idx M
else
    % openadapt coil sensitivities
    [~, ~, cmaps_open] = openadapt_jih(images.pc1, 8, 8, 1);
    cmaps              = conj(cmaps_open);
end

%% Density compensation
% For an explanation of why we do this, see this paper by Thomas
% Benkert: "Optimization and validation of accelerated golden-angle 
% radial sparse MRI reconstruction with self-calibrating GRAPPA 
% operator gridding". MRM 2018.

DATA = DATA .* permute(repmat(sqrt(squeeze(dcf)), [1 1 NCoils]), [1 3 2]);

%% Dot product matching reconstruction
if params_reco.DirectMatching

    % reco undersampled images
    temp_cmap   = conj(cmaps(center_pix,center_pix,:));
    temp_direct = complex(zeros(Nxy, Nxy, NR, 'single'));
    parfor j=1:NR
        FTu = NUFFT(ktraj(:,:,j), dcf(:,:,j), [0 0], [Nxy Nxy]);
        temp_direct(:,:,j) = sum(temp_cmap.*(FTu' * reshape(DATA(:,:,j), [NRead 1 NCoils])), 3);
    end
    clear temp_cmap;

    % zero interpolation filling
    if params_reco.zero_params.on_off
        temp_direct = permute(mg_zero_filling(permute(temp_direct, [3,1,2]), params_reco.zero_params), [2,3,1]);
    end

    % save directly reconstructed images
    if ~isfield(params_reco, 'save_direct')
        params_reco.save_direct = false;
    end
    if params_reco.save_direct
        images.direct = temp_direct;
    end

    % sliding window reconstruction -> useful for analyzing motion
    if isfield(params_reco, 'NSeg')
        NRperSeg              = NR / params_reco.NSeg;
        images.sliding_window = zeros(size(temp_direct,1), size(temp_direct,2), params_reco.NSeg, 'single');
        for j=1:params_reco.NSeg
            tt = (j-1)*NRperSeg+1 : j*NRperSeg;
            images.sliding_window(:,:,j) = mean(temp_direct(:,:,tt),3);
            clear tt;
        end
    end    

    % dictionary matching
    match.direct = mrf_patternmatch(temp_direct, dict, params_dict.look_up, params_dict.look_up_names, [], params_reco.NBlocks);
    clear temp_direct;

end

%% Dot product matching using compressed dictionary and subspace images
if params_reco.DirectMatching_SVD
    images.SVD = lowrankMRF2D_adjoint(DATA, cmaps, proj_id, dict_phi, FT, Nid);
    images.SVD = images.SVD(center_pix,center_pix,:);
    if params_reco.zero_params.on_off
        images.SVD = permute(mg_zero_filling(permute(images.SVD, [3,1,2]), params_reco.zero_params), [2,3,1]);
    end
    match.SVD = mrf_patternmatch(images.SVD, dict_svd, params_dict.look_up, params_dict.look_up_names, [], params_reco.NBlocks);
end

%% Low-Rank Reconstruction
% See these papers for more details...
% Gastao Cruz, MRM 2019. "Sparsity and locally low rank regularization for MR fingerprinting".
% Jesse Hamilton, NMR Biomed 2019. "Simultaneous multislice cardiac magnetic resonance fingerprinting using low rank reconstruction".
if params_reco.LowRank
    Et      = @(x)lowrankMRF2D_adjoint(x, cmaps, proj_id, dict_phi, FT, Nid);          % adjoint operator (spiral k-space to image domain)
    E       = @(x)lowrankMRF2D_forward(x, dict_phi, FT, proj_id, cmaps, [NRead Nid]);  % forward operator (image domain to spiral k-space)
    y0      = E(Et(DATA));
    unitv   = sum(abs(DATA(:))) / sum(abs(y0(:)));
    x0      = Et(DATA);
    scaling = max(abs(x0(:)));
    x0      = x0/scaling;
    DATA    = DATA/scaling;

    % start low rank reconstruction
    params_LR.t0 = unitv;
    images.LR    = nonlinearCGDescent(x0*0, [], E, Et, DATA, params_LR);
    images.LR    = images.LR(center_pix, center_pix,:);
    if params_reco.zero_params.on_off
        images.LR = permute(mg_zero_filling(permute(images.LR, [3,1,2]), params_reco.zero_params), [2,3,1]);
    end
    match.LR = mrf_patternmatch(images.LR, dict_svd, params_dict.look_up, params_dict.look_up_names, [], params_reco.NBlocks);

    % save uncompressed images
    Nxy = size(images.LR, 1);
    if params_LR.save_uncomp
        temp_uncomp = zeros(Nxy, Nxy, NR, 'single');
        parfor y=1:Nxy
        for    x=1:Nxy
            temp_pc = repmat(squeeze(images.LR(y,x,:)), 1, NR);
            temp_s  = sum(dict_phi.' .* temp_pc);
            temp_uncomp(y,x,:) = temp_s(:);
        end
        end
        images.LR_uncomp = temp_uncomp;
    end
    clear FT Et E unitv y0 x0 scaling temp_uncomp;
end

%% permute images
fnames = fieldnames(images);
for j=1:numel(fnames)
    eval(['images.' fnames{j} ' = permute(images.' fnames{j} ', [3,1,2]);']);
end

end
