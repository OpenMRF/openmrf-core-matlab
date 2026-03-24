function [ T12p_Map, S0_Map, R2_Map, RSS_Map, mask_fit] = mg_map_T12p(Images, tSL, mask_fit, mod_fit)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 23.03.2026

% ----- Input: -----
% Images:   N x Ny x Nx, reconstructed images
% tSL:      [s] spin-lock times or echo times
% mask_fit: Ny x Nx, fit mask
% mod_fit:  0: toolbox, 1 fminsearch, 2 fmincon

% ----- Output: -----
% T12p_Map: Ny x Nx, [s] parameter map for T1p, T2p or T2
% S0_Map:   Ny x Nx, fit result for S0
% R2_Map:   Ny x Nx, R2 fit results
% RSS_Map:  Ny x Nx, residula sum of squres
% mask_fit: Ny x Nx, fit mask

% ----- fit model: -----
% S(tSL) = S0 * exp(-tSL/T1p)

%% rotate images to real axis
if ~isreal(Images(:))
    Images = real(Images .* exp(-1i*repmat(angle(Images(1,:,:)),[size(Images,1),1,1]))); % use the first image as the reference
end

%% default parameters
[NtSL, Ny, Nx] = size(Images);

if isempty(mask_fit)
    mask_fit = mg_get_mask_fit(squeeze(abs(mean(Images))), 'holes');
end
if nargin<4
    mod_fit = [];
end
if isempty(mod_fit)
    mod_fit = 1; %  0: toolbox   1: fminsearch   2: fmincon
end

%% save results in arrays
T12p_Map = zeros(Ny, Nx);
S0_Map   = zeros(Ny, Nx);
R2_Map   = zeros(Ny, Nx);
RSS_Map  = zeros(Ny, Nx);

%% ----- prefit: find start values for S0 and T1p -----
ydata_pre = zeros(NtSL,1);
for j=1:NtSL
    temp1(:,:)     = Images(j,:,:);
    temp2(:)       = temp1(mask_fit==1);
    ydata_pre(j,1) = mean(temp2);
end
clear temp1 temp2;

temp1(:,:)            = Images(1,:,:);
temp2(:)              = temp1(mask_fit==1);
S0_start              = mean(temp2);
T1p_start             = tSL(end)*0.75;
[S0_start, T1p_start] = mg_fit_exp(tSL, ydata_pre, S0_start, T1p_start, mod_fit);
clear temp1 temp2
        
%% pixel:pixel fit (parallel)
parfor y=1:Ny
for x=1:Nx
    if (mask_fit(y,x) == 1)
        xdata               = tSL;
        ydata               = Images(:,y,x);
        [S0, T12p, R2, RSS] = mg_fit_exp(xdata, ydata, S0_start, T1p_start, mod_fit);
        S0_Map(y,x)         = S0;
        T12p_Map(y,x)       = T12p;
        R2_Map(y,x)         = R2;
        RSS_Map(y,x)        = RSS;
    end
end
end

end

