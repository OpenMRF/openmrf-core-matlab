% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% case: quadratic, isotropic
if ( ~isfield(FOV, 'Nx') && isfield(FOV, 'Nxy') )
    FOV.Nx = FOV.Nxy;
end
if ( ~isfield(FOV, 'Ny') && isfield(FOV, 'Nxy') )
    FOV.Ny = FOV.Nxy;
end
if ( ~isfield(FOV, 'fov_x') && isfield(FOV, 'fov_xy') )
    FOV.fov_x = FOV.fov_xy;
end
if ( ~isfield(FOV, 'fov_y') && isfield(FOV, 'fov_xy') )
    FOV.fov_y = FOV.fov_xy;
end

%% spatial resolution
if ~isfield(FOV, 'dx') 
    FOV.dx = FOV.fov_x / FOV.Nx;
end
if ~isfield(FOV, 'dy') 
    FOV.dy = FOV.fov_y / FOV.Nx;
end

%% define FOV and encoding steps in z direction: only important for "stacked" 3D readouts
if ~isfield(FOV, 'fov_z') 
    FOV.fov_z = FOV.dz;
end
if ~isfield(FOV, 'Nz') 
    FOV.Nz = 1;
end

%% shift slice position: default 0
if ~isfield(FOV, 'z_offset') 
    FOV.z_offset = 0; % positive/negative -> food/head direction
end

%% reorder struct field names
temp_names = fieldnames(FOV);
[~, idx]   = sort(temp_names);
FOV        = orderfields(FOV, temp_names(idx));
clear temp_names idx;

%% standard orientation
% this might change if you rotate at the scanner!
FOV.orientation     = 'axial';
FOV.read_direction  = 'x';
FOV.phase_direction = 'y';
FOV.Rot_Mat         = [1 0 0; 0 1 0; 0 0 1];