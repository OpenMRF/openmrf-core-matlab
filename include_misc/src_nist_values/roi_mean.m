function mu = roi_mean(map, xc, yc, r, pv_tresh, r2_map, r2_tresh)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 23.03.2026
    
    % ----- Input: -----
    % map:       Ny x Nx, parameter map from which the ROI mean is calculated
    % xc:        x-coordinate of the circular ROI center [pixel]
    % yc:        y-coordinate of the circular ROI center [pixel]
    % r:         radius of the circular ROI [pixel]
    % pv_tresh:  partial-volume threshold for ROI inclusion (default: 0.25)
    % r2_map:    Ny x Nx, optional R2 map for excluding low-quality pixels
    % r2_tresh:  R2 threshold for ROI masking (default: 0.95)

    % ----- Output: -----
    % mu:        weighted mean value of map within the circular ROI

    % defaults
    if isempty(pv_tresh)
        pv_tresh = 0.25;
    end
    if isempty(r2_tresh)
        r2_tresh = 0.95;
    end

    % define ROI with partial volume weighting
    N        = size(map,1);
    [xg, yg] = meshgrid(1:N, 1:N);    
    s        = 20;
    offs     = ((0:s-1)+0.5)/s - 0.5;
    [ox, oy] = meshgrid(offs, offs);
    ox       = ox(:)';
    oy       = oy(:)';    
    dx       = xg(:) - xc;
    dy       = yg(:) - yc;    
    in       = ((dx + ox).^2 + (dy + oy).^2) <= r^2;
    ROI      = reshape(mean(in, 2), N, N);
    ROI      = ROI / max(ROI(:));

    % calculate weighted mean depending on pv_tresh and r2_tresh
    if isempty(r2_map) % no R2 filtering        
        vals    = map(ROI > pv_tresh);
        weights = ROI(ROI > pv_tresh);

    else % adaptive R2 thresholding to keep at least 5 voxels        
        r2_tresh_init = r2_tresh;
        min_voxels    = 5;
        while true
            mask    = (ROI > pv_tresh) .* (r2_map >= r2_tresh);
            vals    = map(mask==1);
            weights = ROI(mask==1);
            if numel(vals) >= min_voxels || r2_tresh <= 0
                break;
            end
            r2_tresh = r2_tresh - 0.01;
        end
        if r2_tresh < r2_tresh_init
            warning('roi_mean:AdjustedR2Thresh', ...
                    'R2 threshold was reduced from %.2f to %.2f to retain at least %d voxels in the ROI.', ...
                    r2_tresh_init, r2_tresh, min_voxels);
        end
    end

    mu = sum(vals .* weights) / sum(weights);

end