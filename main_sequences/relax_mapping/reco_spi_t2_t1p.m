% reconstruction of pulseq spiral mapping sequences
% use for T2/T1p mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% reconstruction of spiral datasets
[Images, PULSEQ, study_info] = SPI_reco();

% fit mask
Nimages = size(Images,1);
[mask, mask3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', Nimages);

%% T2 or T1p mapping

% rotate images to real axis: use the first image as the reference
Images = real(Images .* exp(-1i*repmat(angle(Images(1,:,:)),[Nimages,1,1])));

if isfield(PULSEQ, 'T2')
    TE = PULSEQ.T2.prep_times;
    params.TE = TE;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, TE, mask );
end
if isfield(PULSEQ, 'MLEV')
    t2prep = PULSEQ.MLEV.t2_prep_times;
    params.t2prep = t2prep;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, t2prep, mask );
end
if isfield(PULSEQ, 'SL')
    tSL = PULSEQ.SL.tSL;
    params.tSL = tSL;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, tSL, mask );
end
if isfield(PULSEQ, 'ADIASL')
    tSL = PULSEQ.ADIASL.adia_prep_times;
    params.tSL = tSL;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, tSL, mask );
end

%% vis results
t2cmp  = get_cmp('T2', 1000, 1);

figure;
tiledlayout(1, 2);

nexttile; 
imagesc(T12p_Map, [0 1]); 
axis image; axis off; 
if isfield(PULSEQ, 'T2') || isfield(PULSEQ, 'MLEV')
    title('T2 [s]');
else
    title('T1\rho [s]')
end
colormap(t2cmp);
colorbar;

nexttile;
imagesc(R2_Map, [0.8 1]); 
axis image; axis off; 
title('R^2');
colorbar;

colormap(turbo(1000));

%%
res.T12p_Map = T12p_Map; res.Images = Images;
save_study_results(res, params, study_info, 'prompt', 1);
