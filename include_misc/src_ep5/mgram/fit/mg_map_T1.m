function [T1_Map, M0_Map, Eff_Map, R2_Map, mask_fit] = mg_map_T1(Images, TI, mask_fit)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 23.03.2026

% ----- Input: -----
% Images:   N x Ny x Nx, reconstructed images
% TI:       [s] inversion recovery times
% mask_fit: Ny x Nx, fit mask

% ----- Output: -----
% T1p_Map:  Ny x Nx, [s] parameter map for T1
% M0_Map:   Ny x Nx, fit result for M0
% Eff_Map:  Ny x Nx, inversion efficiency map
% R2_Map:   Ny x Nx, R2 fit results
% mask_fit: Ny x Nx, fit mask

% ----- fit model: -----
% fit model accoring to doi.org/10.1002/mrm.28779
% S = P(1) * (1-(1+P(2))*exp(-TI/T1))

%% rotate images to real axis
if ~isreal(Images(:))
    Images = real(Images .* exp(-1i*repmat(angle(Images(end,:,:)),[size(Images,1),1,1]))); % use the last image as the reference
end

%% default parameters
if nargin<3
    mask_fit = [];
end
if isempty(mask_fit)
    mask_fit = mg_get_mask_fit(squeeze(abs(mean(Images))), 'holes');
end
if size(TI,1) < size(TI,2)
    TI = TI';
end

%% pixel:pixel fit

Ny = size(Images,2);
Nx = size(Images,3);

T1_Map    = zeros(Ny,Nx);
M0_Map    = zeros(Ny,Nx);
Eff_Map   = zeros(Ny,Nx);
R2_Map    = zeros(Ny,Nx);

parfor y = 1:Ny
   for x = 1:Nx
   if mask_fit(y,x)==1

       data    = squeeze(Images(:,y,x));
       P_start = [abs(data(end)), 1.0, 1.0];

       opts   = optimset('Display','off');       
       FITfun = @(P) P(1) * (1-(1+P(2))*exp(-TI/P(3)));
       RSSfun = @(P) sum((data - FITfun(P)).^2);
       P_fit  = fminsearch(RSSfun, P_start, opts);

       RSS    = RSSfun(P_fit);
       TSS    = sum( (data-mean(data)).^2 );
       R2     = 1 - RSS/TSS;

       M0_Map(y,x)  = P_fit(1); 
       Eff_Map(y,x) = P_fit(2);
       T1_Map(y,x)  = P_fit(3);
       R2_Map(y,x)  = R2;
       
   end
   end
end

end

