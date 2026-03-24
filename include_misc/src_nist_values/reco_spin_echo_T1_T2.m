% reconstruction of T1 and T2 weighted spin-echo images
% use for T1/T2 mapping of the NIST phantom
% mgram; V1; 19.02.2026

%% start reco
clear

study = 1;

switch study

    case 1  % rawdata: Avanto 1.5T (Wuerzburg)
    vendor         =  'Siemens';    
    ref_name       =  '1.5T_MnCl2_Wuerzburg';
    study_path     = 'Q:/data/Pulseq/Rawdata/tomgr/Avanto/250912_OpenMRF_NIST';
    study_names_t1 = { 'meas_MID00845_FID16202_se_ir_35.dat',
                       'meas_MID00846_FID16203_se_ir_60.dat',
                       'meas_MID00847_FID16204_se_ir_95.dat',
                       'meas_MID00848_FID16205_se_ir_160.dat',
                       'meas_MID00849_FID16206_se_ir_260.dat',
                       'meas_MID00850_FID16207_se_ir_430.dat',
                       'meas_MID00851_FID16208_se_ir_710.dat',
                       'meas_MID00852_FID16209_se_ir_1180.dat',
                       'meas_MID00853_FID16210_se_ir_1940.dat',
                       'meas_MID00854_FID16211_se_ir_3200.dat'};
    study_names_t2 = { 'meas_MID00855_FID16212_se_te_10.dat',
                       'meas_MID00856_FID16213_se_te_15.dat',
                       'meas_MID00857_FID16214_se_te_20.dat',
                       'meas_MID00858_FID16215_se_te_30.dat',
                       'meas_MID00859_FID16216_se_te_40.dat',
                       'meas_MID00860_FID16217_se_te_60.dat',
                       'meas_MID00861_FID16218_se_te_85.dat',
                       'meas_MID00862_FID16219_se_te_120.dat',
                       'meas_MID00863_FID16220_se_te_175.dat',
                       'meas_MID00864_FID16221_se_te_250.dat'};

    case 2  % rawdata: United Imaging, uMR790, 1.5T (Houston)
    vendor      =  'UI'; 
    ref_name    =  '1.5T_MnCl2_Houston';
    study_path  = 'Q:/data/Pulseq/Rawdata/tomgr/United_Imaging_1,5T/260218_OpenMRF_NIST/mat/';
    study_name  = 'nist_ref_rawdata.mat';

end

%% load rawdata

switch vendor
    case 'Siemens' % load rawdata from twix object
        for j=1:numel(study_names_t1)
            twix_obj = mapVBVD(fullfile(study_path, study_names_t1{j}), 'ignoreSeg', 'removeOS');
            twix_obj = twix_obj{end};
            twix_obj = squeeze(twix_obj.image());
            rawdata_T1(j,:,:,:) = permute(twix_obj,[2,3,1]); % TI x coils x Ny x Nx
            clear twix_obj;
        end
        for j=1:numel(study_names_t2)
            twix_obj = mapVBVD(fullfile(study_path, study_names_t2{j}), 'ignoreSeg', 'removeOS');
            twix_obj = twix_obj{end};
            twix_obj = squeeze(twix_obj.image());
            rawdata_T2(j,:,:,:) = permute(twix_obj,[2,3,1]); % TE x coils x Ny x Nx
            clear twix_obj;
        end

    case 'UI' % load rawdata from mat files        
        load([study_path study_name]);
        rawdata_T1 = permute(irse_rawdata, [1 4 2 3]); % TI x coils x Ny x Nx
        rawdata_T2 = permute(sese_rawdata, [1 4 2 3]); % TE x coils x Ny x Nx
        rawdata_T1 = (rawdata_T1(:,:,:,1:2:end) + rawdata_T1(:,:,:,2:2:end)) / 2;
        rawdata_T2 = (rawdata_T2(:,:,:,1:2:end) + rawdata_T2(:,:,:,2:2:end)) / 2;
        clear irse_rawdata sese_rawdata;
end

%% reconstruct images: ifft
Images_coils_T1 = kspace2image(rawdata_T1);
Images_coils_T2 = kspace2image(rawdata_T2);
clear rawdata_T1 rawdata_T2;

% calculate cmaps via espirit
[cmaps_T1, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils_T1(end,:,:,:)), [], [], [], []);  % use the first image as the reference
[cmaps_T2, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils_T2(1,:,:,:)),   [], [], [], []);  % use the last image as the reference

% calculate coil combined images
for j=1:size(Images_coils_T1, 1)
    Images_T1(j,:,:,:) = squeeze(sum(squeeze(Images_coils_T1(j,:,:,:)) .* conj(cmaps_T1)));
end
for j=1:size(Images_coils_T2, 1)
    Images_T2(j,:,:,:) = squeeze(sum(squeeze(Images_coils_T2(j,:,:,:)) .* conj(cmaps_T2)));
end
clear j Images_coils_T1 Images_coils_T2 cmaps_T1 cmaps_T2;

% rotate images to real axis
Images_T1 = real(Images_T1 .* exp(-1i*repmat(angle(Images_T1(end,:,:)),[size(Images_T1,1),1,1]))); % use the last image as the reference
Images_T2 = real(Images_T2 .* exp(-1i*repmat(angle(Images_T2(1,:,:)),[size(Images_T2,1),1,1])));   % use the first image as the reference

% zero interpolation filling
zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;
Images_T1 = mg_zero_filling(Images_T1, zero_params);
Images_T2 = mg_zero_filling(Images_T2, zero_params);

%% T1 and T2 mapping

mask_fit = mg_get_mask_fit(squeeze(mean(abs(Images_T1)))+squeeze(mean(abs(Images_T2))), 'holes');

TI = [35 60 95 160 260 430 710 1180 1940 3200] *1e-3;
[T1_Map, ~, ~, R2_Map_T1] = mg_map_T1( double(real(Images_T1)), TI, mask_fit);

TE = [10 15 20 30 40 60 85 120 175 250] *1e-3;
[T2_Map, ~, R2_Map_T2] = mg_map_T12p( double(real(Images_T2)), TE, mask_fit );

%% select roi centers
N = 14; % number of inserts in NIST phantom

t1cmp  = get_cmp('T1', 1000, 1);
t2cmp  = get_cmp('T2', 1000, 1);

[xlims, ylims] = meshgrid(1:size(mask_fit,1), 1:size(mask_fit,2));
xlims          = xlims(mask_fit==1);
ylims          = ylims(mask_fit==1);
xlims          = [min(xlims), max(xlims)] + 3*[-1 1];
ylims          = [min(ylims), max(ylims)] + 3*[-1 1];

[xc, yc] = pick_phantom_points( N, ...
                                T1_Map, T1_Map, T1_Map, T1_Map, ...
                                [0 0.5], [0 4], [0 0.2], [0 3], ...
                                t1cmp, t1cmp, t2cmp, t2cmp, ...
                                xlims, ylims );

%% calc mean in ROIs
r  = 3;
for j = 1:N
    T1(j,1) = roi_mean(T1_Map, xc(j), yc(j), r, 0.25, R2_Map_T1, 0.95);
    T2(j,1) = roi_mean(T2_Map, xc(j), yc(j), r, 0.25, R2_Map_T2, 0.95);
end

%% vis results
figure;
ax1 = subplot(2,2,1);
imagesc(T1_Map, [0 0.4]); axis image; axis off; colormap(ax1, t1cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);
ax2 = subplot(2,2,2);
imagesc(T1_Map, [0 3]); axis image; axis off; colormap(ax2, t1cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);
ax3 = subplot(2,2,3);
imagesc(T2_Map, [0 0.05]); axis image; axis off; colormap(ax3, t2cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);
ax4 = subplot(2,2,4);
imagesc(T2_Map, [0 2]); axis image; axis off; colormap(ax4, t2cmp); hold on; plot(xc, yc, 'gx'); plot(xc, yc, 'r.');
xlim(xlims); ylim(ylims);

%% compare to nominal values
[T1_nom, T2_nom] = NIST_references(ref_name(1:10));

figure()
subplot(1,2,1)
loglog(T1_nom, T2_nom, 'go', 'MarkerSize', 10)
hold on
loglog(T1, T2, 'k.', 'MarkerSize', 20)
xlabel('T1 [s]')
ylabel('T2 [s]')
legend('nom', 'meas', 'Location', 'best')
axis square
grid on
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'LineWidth', 2, 'Fontsize', 16); 

subplot(1,2,2)
semilogx(T1_nom, (T1-T1_nom)./T1_nom*100, 'r.', 'MarkerSize', 20)
hold on
semilogx(T2_nom, (T2-T2_nom)./T2_nom*100, 'b.', 'MarkerSize', 20)
yline(0, 'k--')
xlabel('T1/T2 [s]')
ylabel('deviation [%]')
legend('T1', 'T2', 'Location', 'best')
ylim([-50 50])
axis square
grid on
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'LineWidth', 2, 'Fontsize', 16); 

%% save results
save([ref_name '.mat'], 'T1_Map', 'T2_Map', 'R2_Map_T1', 'R2_Map_T2', 'T1', 'T2', 'xc', 'yc', 'r');