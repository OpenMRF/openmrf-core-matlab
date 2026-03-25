%%
clear
pulseq_init();

% params: field of view
FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 48;          % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 96   *1e-3;  % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

% params: basic
SPI.mode_2D_3D = '3D_stacked';     % '2D', '3D' or '3D_stacked'
SPI.TR         = 20 *1e-3; % [s]  repetition time, 0 for minimization
SPI.TE         = 0 *1e-3;  % [s]  echo time, 0 for minimization
SPI.NR         = 100;      % [ ]  number of repetitions
SPI.Ndummy     = 0;        % [ ]  initial dummy loops
SPI.adcBW      = 100 *1e3; % [Hz] desired receiver bandwidth

% params: slice excitation
SPI.exc_mode      = 'sinc';     % 'sinc' or 'sigpy_SLR'
SPI.exc_shape     = 'ex';       % only for sigpy: 'st' or 'ex' 
SPI.exc_time      = 1.0 *1e-3;  % [s] excitation time
SPI.exc_tbw       = 6;          % [ ] time bandwidth product
SPI.exc_fa_mode   = 'equal';    % 'equal',  'ramped',  'mrf'  
SPI.exc_flipangle = 20 *pi/180; % [rad] const FA or start of FA ramp -> set [] for auto mode

% params: gradient spoiling & rf spoiling
SPI.spoil_nTwist  = 8;           % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode = 'lin';       % 'lin' or 'quad'
SPI.spoil_rf_inc  = 117 *pi/180; % [rad] rf phase increment

% k-space geometry params
id_test = 1;
switch id_test
case 1
SPI.geo.design        = 'spiral';
SPI.geo.design_fun    = 'hargreaves';
SPI.geo.N_interleaves = 24;
SPI.geo.FOV_coeff     = [1 -0.5];
SPI.geo.Ns            = 1e5;
SPI.geo.ds            = 1e-3;
SPI.geo.flag_rv       = 0;

case 2        
SPI.geo.design   = 'rosette';
SPI.geo.N_lobes  = 17;
SPI.geo.Ns       = 1e5;
SPI.geo.ds       = 1e-3;
SPI.geo.flag_rv  = 0;

case 3        
SPI.geo.design   = 'radial';
SPI.geo.t_adc    = 3.2 *1e-3;

case 4        
SPI.geo.design     = 'radialHalf';
SPI.geo.design_fun = 'minTimeTrap';

end

% k-space projection params
SPI.proj.mode = 'RoundGoldenAngle'; % 'Equal2Pi' 'RoundGoldenAngle'
SPI.proj.Nid  = 48; % number of unique projections

% gradient & slew rate limitation factors
SPI.lim_exc_slew   = 1.0; % slice excitation
SPI.lim_reph_grad  = 1.0; % slice rephasing
SPI.lim_reph_slew  = 1.0; % slice rephasing
SPI.lim_read_grad  = 1.0; % readout
SPI.lim_read_slew  = 1.0; % readout
SPI.lim_spoil_grad = 1.0; % slice spoiling
SPI.lim_spoil_slew = 1.0; % slice spoiling

% calculate SPI pulse objects
[SPI, ktraj_ref] = SPI_init(SPI, FOV, system, 1);

%% add sequence blocks
for loop_kz = 1 : FOV.Nz
    for loop_NR = 1-SPI.Ndummy : 10%SPI.NR
        SPI_add();
    end
end

%% plot sequence
seq.plot();
[timings_ok, error_report] = seq.checkTiming;

%% check full k-space trajectory
ktraj_adc_full = pulseq_get_ktraj(seq, 2);

%% check timings and stimulation
if (timings_ok)
    disp('   Timing check passed successfully');
else
    disp(   'Timing check failed! Error listing follows:');
    fprintf([error_report{:}]);
    fprintf('\n');
end

[pns_ok, pns_norm, pns_comp] = PNS_sim(seq, 'axial');