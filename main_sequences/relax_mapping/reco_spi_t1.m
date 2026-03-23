% reconstruction of pulseq spiral mapping sequences
% use for T1 mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% reconstruction of spiral datasets
[Images, PULSEQ, study_info] = SPI_reco();

% fit mask
NInv = size(Images,1);
[mask, mask3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', NInv);

%% T1 mapping

% rotate images to real axis: use the last image as the reference
Images = real(Images .* exp(-1i*repmat(angle(Images(end,:,:)),[NInv,1,1])));

TI = PULSEQ.INV.inv_rec_time;
[T1_Map, M0_Map, Eff_Map, R2_Map] = mg_map_T1(real(Images), TI, mask);

%% vis results
t1cmp  = get_cmp('T1', 1000, 1);

figure;
tiledlayout(1, 3)

nexttile;
imagesc(T1_Map, [0 2]); 
axis image; axis off; 
title('T1 [s]');
colormap(t1cmp);
colorbar;

nexttile;
imagesc(R2_Map, [0.8 1]); 
axis image; axis off; 
title('R^2')
colorbar;

nexttile;
imagesc(Eff_Map, [0.5, 1]); 
axis image; axis off; 
title('Inversion efficiency')
colorbar;

colormap(turbo(1000));

%%
res.T1_Map = T1_Map; res.Images = Images;
params.TI = TI;
save_study_results(res, params, study_info);