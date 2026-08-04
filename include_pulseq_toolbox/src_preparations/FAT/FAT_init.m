function FAT = FAT_init(FAT, FOV, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% FAT saturation
if strcmp(FAT.mode, 'on')

    %  http://mriquestions.com/uploads/3/4/5/7/34572113/fat-water_review_jmri.21895.pdf
    FAT.ppm   = -3.45;        % [ppm] chemical shift of fat fraction
    FAT.df0   = FAT.ppm *1e-6 * system.B0 * system.gamma;
    if ~isfield(FAT, 'tExc')
        FAT.tExc  = 8  *1e-3;     % [s] pulse duration
    end
    if ~isfield(FAT, 'alpha')
        FAT.alpha = 100 *pi/180;  % [deg] flip angle
    end
    if ~isfield(FAT,'rf_type')
        FAT.rf_type = 'slr'; % default
    end
    if strcmp(FAT.rf_type, 'slr')
        if ~isfield(FAT, 'rf_d1')
            FAT.rf_d1 = 0.01;
        end
        if ~isfield(FAT, 'rf_d2')
            FAT.rf_d2 = 0.01;
        end
    end
    if ~isfield(FAT, 'rf_bw_f')
        FAT.rf_bw_f = 1;
    end

    if strcmp(FAT.rf_type, 'slr')
        FAT.rf = SIGPY_SLR(  FAT.alpha, ...
                             FAT.tExc,  ...
                             0,  ...
                             abs(FAT.tExc * FAT.df0) * FAT.rf_bw_f ,  ...
                             'sat',  ...
                             'ls',  ...
                             FAT.rf_d1, FAT.rf_d2, ...
                             0 ,  ...
                             1e6,  ...
                             system);
        FAT.rf.freqOffset = FAT.df0;
    end

    if strcmp(FAT.rf_type, 'gauss')
        FAT.rf    = mr.makeGaussPulse(  FAT.alpha, system, ...
                                        'Duration', FAT.tExc, ...
                                        'bandwidth', abs(FAT.df0) * FAT.rf_bw_f, ...
                                        'freqOffset', FAT.df0, ...
                                        'SliceThickness', 1e6, ...
                                        'PhaseOffset', 0, ...
                                        'use', 'saturation' );
    end
    
    if strcmp(FAT.rf_type, 'sinc')
        FAT.rf    = mr.makeSincPulse(   FAT.alpha, system, ...
                                        'Duration', FAT.tExc, ...
                                        'timeBwProduct', abs(FAT.tExc * FAT.df0) * FAT.rf_bw_f, ...
                                        'freqOffset', FAT.df0, ...
                                        'SliceThickness', 1e6, ...
                                        'PhaseOffset', 0 , ...
                                        'use', 'saturation' );
    end

    % calculate crusher objects
    if ~isfield(FAT, 'crush_nTwists_x')
        FAT.crush_nTwists_x = 2;   % [] number of 2pi twists in x direction
    end
    if ~isfield(FAT, 'crush_nTwists_y')
        FAT.crush_nTwists_y = 2;   % [] number of 2pi twists in y direction
    end
    if ~isfield(FAT, 'crush_nTwists_z')
        FAT.crush_nTwists_z = 5.1;   % [] number of 2pi twists in z direction
    end
    if ~isfield(FAT, 'crush_lim_grad')
        FAT.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
    end    
    if ~isfield(FAT, 'crush_lim_slew')
        FAT.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
    end
    [FAT.gx_crush, FAT.gy_crush, FAT.gz_crush] = CRUSH_x_y_z(FAT.crush_nTwists_x, FAT.crush_nTwists_y, FAT.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, FAT.crush_lim_grad, FAT.crush_lim_slew, system);
    FAT.tcrush = mr.calcDuration(FAT.gx_crush, FAT.gy_crush, FAT.gz_crush);

end

end