function [SPI, ktraj_ref] = SPI_init(SPI, FOV, system, flag_plot)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ----------------- readout: spiral (SPI) -----------------
% ---------------------------------------------------------

% Comment: Originally, this was a pure spiral readout module. However,
% after years of implementing new features, the SPI package is more like
% an universal readout module, which can be used for different acquistions
% e.g. radial, spiral, rosette, cones, seifferts ...

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 10.03.2026

% -------------------------------------------------------------------
% basic checks
% calculate rf pulses and slice selection gradient
% calculate combined slice rephaser and partition gradients
% calculate k-space trajectory and optimized gradients for xyz-readout & xyz-rewinder
% calculate ADC duration and ADC object
% calculate gx gy and gz gradient objects
% calculate ADC padding (adjust ADC delay for radial timing)
% calculate projection angles and rotation events
% calculate slice spoiler or combined rewinder gradients
% set rf spoiling
% calculate Skope field cam objects
% adapt timings for TRs and TEs
% export reference k-space trajectory for reconstruction
% -------------------------------------------------------------------

% ----- INPUT -----
% SPI:       structured object containing initial readout parameters
% FOV:       structured object containing the field of view parmeters
% system:    structured object containing the MRI system information
% flag_plot: plot the output k-space trajectory (optional)

% ----- OUTPUT -----
% SPI:       structured object containing the final readout parameters and sequence objects
% ktraj_ref: [3 x Nadc] reference k-space trajectory for SPI.proj.phi=0 and SPI.theta=0

% optional import parameters for MRF use via SPI.mrf_import:
% phi, dphi, phi_id -> rotation angles of spiral interleaves
% FAs -> flip agnles
% TRs -> repetition times
% TEs -> echo times
% rf_phase -> rf spoiling increments
% gxy -> gradient trajecotry
% kxy -> k-space trajecotry

%% MRF IP notice
if isfield(SPI, 'mrf_import')
    mrf_ip_notice();
end

%% basic checks

if ~isfield(SPI, 'Ndummy')
    SPI.Ndummy = 0;
end

% default gradient limits
if ~isfield(SPI, 'lim_exc_slew')
    SPI.lim_exc_slew = 1.0;
end
if ~isfield(SPI, 'lim_reph_grad')
    SPI.lim_reph_grad = 1.0;
end
if ~isfield(SPI, 'lim_reph_slew')
    SPI.lim_reph_slew = 1.0;
end
if ~isfield(SPI, 'lim_read_grad')
    SPI.lim_read_grad = 1.0;
end
if ~isfield(SPI, 'lim_read_slew')
    SPI.lim_read_slew = 1.0;
end
if ~isfield(SPI, 'lim_spoil_grad')
    SPI.lim_spoil_grad = 1.0;
end
if ~isfield(SPI, 'lim_spoil_slew')
    SPI.lim_spoil_slew = 1.0;
end

% check balanced or unbalanced mode
if abs(SPI.spoil_nTwist) > 0
    SPI.spoil_gz_mode = 'unbalanced';
else
    SPI.spoil_gz_mode = 'balanced';
end
if strcmp(SPI.spoil_gz_mode, 'balanced')
    SPI.spoil_nTwist = 0;
end

% check acquisition mode: '2D', '3D' or '3D_stacked'
if ~isfield(SPI, 'mode_2D_3D')
    error('choose SPI.mode_2D_3D: >2D<, >3D< or >3D_stacked<');
end
if ~(strcmp(SPI.mode_2D_3D, '2D') || strcmp(SPI.mode_2D_3D, '3D') || strcmp(SPI.mode_2D_3D, '3D_stacked'))
    error('choose SPI.mode_2D_3D: >2D<, >3D< or >3D_stacked<');
end
if strcmp(SPI.mode_2D_3D, '3D_stacked')
    if ~isfield(FOV, 'Nz')
        error('FOV.Nz is required in a 3D_stacked sequence!');
    end
end
if strcmp(SPI.mode_2D_3D, '2D') || strcmp(SPI.mode_2D_3D, '3D')
    if ~isfield(FOV, 'Nz')
        FOV.Nz = 1;
    end
    if FOV.Nz ~= 1
        error('FOV.Nz (number of z partitions) is 1 for 2D or 3D readouts');
    end
end
if strcmp(SPI.geo.design, 'cones') || strcmp(SPI.geo.design, 'seiffert')
    if strcmp(SPI.mode_2D_3D, '2D') || strcmp(SPI.mode_2D_3D, '3D_stacked')
        error('cones or seiffert readouts are not comaptible with 2D or 3D_stacked mode');
    end
end

%% calculate rf pulses and slice selection gradient

% calculate or import flip angles
switch SPI.exc_fa_mode
    case 'equal'
        SPI.exc_flipangle = SPI.exc_flipangle * ones(SPI.NR,1);
    case 'ramped'
        SPI.exc_flipangle = SPI_calc_rf_ramp(SPI.NR, SPI.exc_flipangle);
    case 'import'
        SPI.exc_flipangle = SPI.mrf_import.FAs;
end

% create rf objects and gz dummy
for j = 1:SPI.NR
    if strcmp(SPI.exc_mode, 'sinc')
        [SPI.rf(j,1), SPI.gz_exc] = mr.makeSincPulse( SPI.exc_flipangle(j), ...
                                                      system, ...
                                                      'Duration', SPI.exc_time,...
                                                      'SliceThickness', FOV.dz, ...
                                                      'timeBwProduct', SPI.exc_tbw, ...
                                                      'apodization', 0.5, ...
                                                      'PhaseOffset', 0, ...
                                                      'maxSlew', system.maxSlew * SPI.lim_exc_slew, ...
										              'use', 'excitation' );
    elseif strcmp(SPI.exc_mode, 'sigpy_SLR')
        [SPI.rf(j,1), SPI.gz_exc] = SIGPY_SLR( SPI.exc_flipangle(j), ...
                                               SPI.exc_time, ...
                                               0, ...
                                               SPI.exc_tbw, ...
                                               SPI.exc_shape, ...
                                               'ls', 0.001, 0.01, 0 , ...
                                               FOV.dz, ...
							   		           system, ...
                                               system.maxSlew * SPI.lim_exc_slew);
    end
    SPI.rf(j).freqOffset = FOV.z_offset * SPI.gz_exc.amplitude;
end
clear j;

%% calculate combined slice rephaser and partition gradients

if ~isfield(SPI, 'reph_duration')
    SPI.reph_duration = [];
end
if isempty(SPI.reph_duration)
    SPI.reph_duration = 0;
end

% calculate slice rephaser for 2D or 3D
if strcmp(SPI.mode_2D_3D, '2D') || strcmp(SPI.mode_2D_3D, '3D')
    if SPI.reph_duration == 0
        SPI.gz_reph(1)    = mr.makeTrapezoid('z', 'Area', -SPI.gz_exc.area/2, 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system);
        SPI.reph_duration = mr.calcDuration(SPI.gz_reph(1));
    else
        SPI.gz_reph(1) = mr.makeTrapezoid('z', 'Area', -SPI.gz_exc.area/2, 'Duration', SPI.reph_duration, 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system);
    end
    SPI.nz_zero = 1;
end

% calculate slice rephaser for 3D_stacked
if strcmp(SPI.mode_2D_3D, '3D_stacked')
    SPI.deltakz = 1 / FOV.fov_z;
    SPI.kz_area = (-FOV.Nz/2 : 1 : FOV.Nz/2-1)' * SPI.deltakz;
    SPI.nz_zero = find(SPI.kz_area==0);
    if SPI.reph_duration == 0
        for j = 1:FOV.Nz
            temp_dur(j) = mr.calcDuration( mr.makeTrapezoid('z', 'Area', -SPI.gz_exc.area/2 + SPI.kz_area(j), 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system) );
        end
        SPI.reph_duration = max(temp_dur);
    end
    for j = 1:FOV.Nz
        SPI.gz_reph(j,1) = mr.makeTrapezoid('z', 'Area', -SPI.gz_exc.area/2 + SPI.kz_area(j), 'Duration', SPI.reph_duration, 'maxGrad', system.maxGrad*SPI.lim_reph_grad, 'maxSlew', system.maxSlew*SPI.lim_reph_slew, 'system', system);
    end
    clear temp_dur;
end
clear j;

%% calculate k-space trajectory and optimized gradients for xyz-readout & xyz-rewinder
SPI.geo.deltak    = 1 / FOV.fov_xy;
SPI.geo.kmax      = FOV.Nxy / 2 * SPI.geo.deltak;
SPI.geo.lim_grad  = SPI.lim_read_grad;
SPI.geo.lim_slew  = SPI.lim_read_slew;
if isfield(SPI.geo, 'FOV_coeff')
    SPI.geo.FOV_coeff = SPI.geo.FOV_coeff * FOV.fov_xy;
end
if ~isfield(SPI.geo, 'flag_plot')
    SPI.geo.flag_plot = 0;
end
switch SPI.geo.design
    case 'spiral'
        SPI.geo = SPI_minTimeSpiral(SPI.geo, system, SPI.geo.flag_plot);
    case 'rosette'
        SPI.geo = SPI_minTimeRosette(SPI.geo, system, SPI.geo.flag_plot);
    case 'radial'
        SPI.geo = SPI_minTimeRadial(SPI.geo, system, SPI.geo.flag_plot);
    case 'radialHalf'
        SPI.geo = SPI_minTimeRadialHalf(SPI.geo, system, SPI.geo.flag_plot);
    case 'cones'
        SPI.geo = SPI_minTimeCones(SPI.geo, system, SPI.geo.flag_plot);
    case 'seiffert'
        SPI.geo = SPI_minTimeSeiffert(SPI.geo, system, SPI.geo.flag_plot);
    case 'import'
        SPI.geo = SPI_importTraj(SPI.geo, system, SPI.geo.flag_plot);
    otherwise
        error('SPI.geo.design is unknown');
end

%% calculate ADC duration and ADC object

% first adc samples might be corrupted by digital filters -> add gradient buffer delay and use adc padding
if ~isfield(SPI, 'adcBufferTime1')
    SPI.adcBufferTime1 = 50 *1e-6; % 50us @500kHz -> 25 samples;  50us @100kHz -> 5 samples
end
% start read gradients with an additional small delay to ensure sampling of the k-space center
if ~isfield(SPI, 'adcBufferTime2')
    SPI.adcBufferTime2 = 10 *1e-6; % 10us @500kHz -> 5 samples;  10us @100kHz -> 1 sample
end
SPI.adcBufferTime = ceil( (SPI.adcBufferTime1+SPI.adcBufferTime2) / system.gradRasterTime ) * system.gradRasterTime;
temp_n_buffer     = ceil( SPI.adcBufferTime / system.gradRasterTime );
SPI.geo.g         = [zeros(temp_n_buffer,3); SPI.geo.g; zeros(1,3)];

% calc adc duration from k-space geometry: search for last kmax visit
temp_kxyz   = sqrt(sum(cumsum(SPI.geo.g).^2,2)) * system.gradRasterTime;
temp_idx    = find(islocalmax(temp_kxyz)==1);
temp_idx    = temp_idx(end);
SPI.adcTime = (temp_idx + temp_n_buffer) * system.gradRasterTime;
clear temp_n_buffer temp_kxyz temp_idx

% special case for spiral trajectories with sampling of the rewinder
if isfield(SPI.geo, 'sample_rewinder') && SPI.geo.sample_rewinder
    SPI.adcTime = size(SPI.geo.g,1) * system.gradRasterTime;
else
    SPI.geo.sample_rewinder = false;
end

% calc ADC object
[SPI.adc, SPI.adcTime, SPI.adcBW, SPI.adcNSamples, SPI.adcDwell] = SPI_calc_adc(SPI.adcTime, SPI.adcBW, system);

%% calculate gx gy and gz gradient objects

% zero padding due to non-zero adc delay
SPI.geo.g = [zeros(round(SPI.adc.delay/system.gradRasterTime),3); SPI.geo.g];

% check adc duration compared to gradient duration
temp_dt = mr.calcDuration(SPI.adc) - size(SPI.geo.g,1)*system.gradRasterTime;
if temp_dt >= 0
    temp_n    = max([ceil(temp_dt/system.gradRasterTime), 1]);
    SPI.geo.g = [SPI.geo.g; zeros(temp_n,3)];
end
clear temp_dt temp_n;

% calc gradient objects with arbitrary waveforms
SPI.gx = mr.makeArbitraryGrad('x', SPI.geo.g(:,1), 'first', 0, 'last', 0, 'system', system);
SPI.gy = mr.makeArbitraryGrad('y', SPI.geo.g(:,2), 'first', 0, 'last', 0, 'system', system);
SPI.gz = mr.makeArbitraryGrad('z', SPI.geo.g(:,3), 'first', 0, 'last', 0, 'system', system);

% update g(t), k(t) and s(t) in geometry struct
SPI.geo.g = [SPI.gx.waveform, SPI.gy.waveform, SPI.gz.waveform];
SPI.geo.k = cumsum(SPI.geo.g) * system.gradRasterTime;
SPI.geo.s = diff(SPI.geo.g)   / system.gradRasterTime;

%% calculate ADC padding (adjust ADC delay for radial timing)

% SPI.adcNPad is a 2x1 vector which defines the adc samples which should be
% used for reconstruction -> SPI.adcNPad(1) ... SPI.adcNPad(2)

% all center-out geometries
warning('OFF', 'mr:restoreShape');
if ~strcmp(SPI.geo.design, 'radial')
    temp_seq = mr.Sequence(system);
    temp_seq.addBlock(SPI.rf(1));
    temp_seq.addBlock(SPI.adc, SPI.gx, SPI.gy, SPI.gz);    
    temp_kxyz      = temp_seq.calculateKspacePP()';
    temp_kxyz      = sqrt(sum(temp_kxyz.^2,2));
    temp_idx       = find(islocalmax(temp_kxyz)==1);
    SPI.adcNPad(1) = round( SPI.adcBufferTime1 / SPI.adc.dwell );
    SPI.adcNPad(2) = temp_idx(end);
    clear temp_seq temp_kxyz temp_idx;

% radial -> let's try to 'hit' the k-space center  
else
    temp_dt = system.adcRasterTime : system.adcRasterTime : 2*SPI.adcDwell;
    for j = 1:numel(temp_dt)
        temp_adc = SPI.adc;
        temp_adc.delay = temp_dt(j);
        temp_seq = mr.Sequence(system);
        temp_seq.addBlock(SPI.rf(1));
        temp_seq.addBlock(temp_adc, SPI.gx, SPI.gy, SPI.gz);        
        temp_kxyz = temp_seq.calculateKspacePP()';
        temp_kxyz = sqrt(sum(temp_kxyz.^2,2));
        temp_idx  = find(islocalmin(temp_kxyz)==1);
        temp_0(j) = temp_kxyz(temp_idx);
    end
    temp_idx = find(temp_0==min(temp_0));
    temp_dt  = temp_dt(temp_idx(1));
    if temp_dt > SPI.adc.dwell
        temp_dt = temp_dt - SPI.adc.dwell;
    end
    SPI.adc.delay = temp_dt;
    temp_seq      = mr.Sequence(system);
    temp_seq.addBlock(SPI.rf(1));
    temp_seq.addBlock(SPI.adc, SPI.gx, SPI.gy, SPI.gz);
    temp_kxyz = temp_seq.calculateKspacePP()';
    temp_kxyz = sqrt(sum(temp_kxyz.^2,2));
    temp_idx  = find(islocalmax(temp_kxyz)==1);
    temp_1    = temp_idx(1);
    temp_2    = temp_idx(2);
    while(temp_kxyz(temp_1)>SPI.geo.kmax)
        temp_1 = temp_1 + 1;
    end
    while(temp_kxyz(temp_2)>SPI.geo.kmax)
        temp_2 = temp_2 - 1;
    end
    SPI.adcNPad(1) = temp_1-1;
    SPI.adcNPad(2) = temp_2+1;
    clear j temp_0 temp_1 temp_2 temp_adc temp_dt temp_idx temp_kxyz temp_seq; 
end

% special case for spiral trajectories with sampling of the rewinder
if SPI.geo.sample_rewinder
    SPI.adcNPad(1) = round( SPI.adcBufferTime1 / SPI.adc.dwell );
    SPI.adcNPad(2) = SPI.adc.numSamples;
end

%% calculate projection angles and rotation events

% calculate or import rotation angles
switch SPI.proj.mode
    case 'Equal2Pi'
        SPI.proj.Nid  = SPI.NR;
        SPI.proj.id   = 1 : SPI.proj.Nid;
        SPI.proj.dphi = 2*pi / SPI.proj.Nid;
        SPI.proj.phi  = SPI.proj.id * SPI.proj.dphi;
        SPI.proj.phi  = wrapTo2Pi(SPI.proj.phi - SPI.proj.phi(1));
    case 'Random2Pi'
        rng("default");
        SPI.proj.Nid     = SPI.NR;
        [~, SPI.proj.id] = sort(rand(1,SPI.proj.Nid));
        SPI.proj.dphi    = 2*pi / SPI.proj.Nid;
        SPI.proj.phi     = SPI.proj.id * SPI.proj.dphi;
        SPI.proj.phi     = wrapTo2Pi(SPI.proj.phi - SPI.proj.phi(1));
    case 'GoldenAngle'
        SPI.proj.Nid  = SPI.NR;
        SPI.proj.id   = 1 : SPI.NR;
        SPI.proj.dphi = 4*pi / (1+sqrt(5));
        SPI.proj.phi  = SPI.proj.id * SPI.proj.dphi;
        SPI.proj.phi  = wrapTo2Pi(SPI.proj.phi - SPI.proj.phi(1));
    case 'RandomGoldenAngle'
        rng("default");
        SPI.proj.Nid     = SPI.NR;
        [~, SPI.proj.id] = sort(rand(1,SPI.proj.Nid));
        SPI.proj.dphi    = 4*pi / (1+sqrt(5));
        SPI.proj.phi     = SPI.proj.id * SPI.proj.dphi;
        SPI.proj.phi     = wrapTo2Pi(SPI.proj.phi - SPI.proj.phi(1));
    case 'RoundGoldenAngle'
        SPI.proj.dphi = 4*pi / (1+sqrt(5));
        SPI.proj.phi  = wrapTo2Pi(SPI.proj.dphi * (1:SPI.NR));
        SPI.proj.dphi = 2*pi / SPI.proj.Nid;
        SPI.proj.phi  = round(SPI.proj.phi/SPI.proj.dphi)' * SPI.proj.dphi;
        SPI.proj.id   = round(SPI.proj.phi/SPI.proj.dphi);
        SPI.proj.id(SPI.proj.id==0) = SPI.proj.Nid;
        SPI.proj.phi  = SPI.proj.id * SPI.proj.dphi;
        SPI.proj.phi  = wrapTo2Pi(SPI.proj.phi - SPI.proj.phi(1));
    case 'FibonacciSphere'
        rng("default");
        SPI.proj.Nid     = SPI.NR;
        [~, SPI.proj.id] = sort(rand(1,SPI.proj.Nid));
        temp_phi         = SPI_fibonacci_sphere(SPI.proj.Nid);
        SPI.proj.phi     = temp_phi(:,1) + 1i * temp_phi(:,2);
        clear temp_phi;
    case 'import'
        SPI.proj.Nid = numel(unique(SPI.mrf_import.phi));  
        SPI.proj.id  = SPI.mrf_import.phi_id;
        SPI.proj.phi = SPI.mrf_import.phi;        
end
SPI.proj.id  = SPI.proj.id(:);  % force column vectors
SPI.proj.phi = SPI.proj.phi(:); % force column vectors

% built unique phi array and check consitency
SPI.proj.phi_unique = unique([SPI.proj.id, SPI.proj.phi], 'rows');
SPI.proj.phi_unique(:,1) = [];
if SPI.proj.Nid ~= numel(SPI.proj.phi_unique)
    error('inconsistent number of unique projections');
end

% calculate roation objects
% note: SPI.proj.phi and SPI.proj.phi_unique can be real or complex valued
% case 1 (real): for 2D readouts, rotations are around z
% case 2 (complex): for 3D readouts, there are two spherical rotation angles
%                   real(phi) -> azimuthal angle
%                   imag(phi) -> polar angle
for j = 1:SPI.proj.Nid
    SPI.rot(j,1) = mr.makeRotation(real(SPI.proj.phi_unique(j)), imag(SPI.proj.phi_unique(j)));
end

%% calculate slice spoiler or combined rewinder gradients

if ~isfield(SPI, 'spoil_duration')
    SPI.spoil_duration = [];
end
if isempty(SPI.spoil_duration)
    SPI.spoil_duration = 0;
end

% final kz value at end of TR
SPI.kz_final = SPI.spoil_nTwist / FOV.fov_z;

% calculate spoiler/rewinder for 2D or 3D
if strcmp(SPI.mode_2D_3D, '2D') || strcmp(SPI.mode_2D_3D, '3D')
    if SPI.spoil_duration==0
        SPI.gz_spoil(1)    = mr.makeTrapezoid('z', 'Area', SPI.kz_final - SPI.gz_exc.area/2, 'maxGrad', system.maxGrad*SPI.lim_spoil_grad, 'maxSlew', system.maxSlew*SPI.lim_spoil_slew, 'system', system);
        SPI.spoil_duration = mr.calcDuration(SPI.gz_spoil(1));
    else
        SPI.gz_spoil(1) = mr.makeTrapezoid('z', 'Area', SPI.kz_final - SPI.gz_exc.area/2, 'Duration', SPI.spoil_duration, 'maxGrad', system.maxGrad*SPI.lim_spoil_grad, 'maxSlew', system.maxSlew*SPI.lim_spoil_slew, 'system', system);  
    end
end

% calculate combined spoiler and rewinder for stacked 3D
if strcmp(SPI.mode_2D_3D, '3D_stacked')
    if SPI.spoil_duration==0
        for j = 1:FOV.Nz
            temp_duration(j) = mr.calcDuration(mr.makeTrapezoid('z', 'Area', SPI.kz_final - SPI.gz_exc.area/2 - SPI.kz_area(j), 'maxGrad', system.maxGrad*SPI.lim_spoil_grad, 'maxSlew', system.maxSlew*SPI.lim_spoil_slew, 'system', system));
        end
        SPI.spoil_duration = max(temp_duration);
    end
    for j = 1:FOV.Nz
        SPI.gz_spoil(j,1) = mr.makeTrapezoid('z', 'Area', SPI.kz_final - SPI.gz_exc.area/2 - SPI.kz_area(j), 'Duration', SPI.spoil_duration, 'maxGrad', system.maxGrad*SPI.lim_spoil_grad, 'maxSlew', system.maxSlew*SPI.lim_spoil_slew, 'system', system);
    end
end
clear temp_duration j;

%% set rf spoiling
switch SPI.spoil_rf_mode
    case 'quad'
        SPI.spoil_rf_pow = 2;
    case 'lin'
        SPI.spoil_rf_pow = 1;
    case 'import'
        SPI.spoil_rf_pow = []; 
    otherwise
        error('unknown rf spoiling mode!')
end

%% calculate Skope field cam objects
if isfield(SPI, 'skope') && isfield(SPI.skope, 'n_seg') && SPI.skope.flag_onoff
    SPI.skope.TR_list = reshape([1:SPI.NR nan(1,ceil(SPI.NR/SPI.skope.n_seg)*SPI.skope.n_seg-SPI.NR)], SPI.skope.n_seg, []).';
    if any(SPI.skope.TR_list(:)<1) || any(SPI.skope.TR_list(:)>SPI.NR)
        error('SPI.skope.TR_list is incompatible with number of spirals!');
    end
    if ~isfield(SPI.skope, 'ttl_dur')
        SPI.skope.ttl_dur = 10e-6;
    end
    if ~isfield(SPI.skope, 'grad_free_dur')
        SPI.skope.grad_free_dur = 300e-6;
    end
    SPI.skope.ttl_dur       = ceil(SPI.skope.ttl_dur/system.blockDurationRaster)       * system.blockDurationRaster;
    SPI.skope.grad_free_dur = ceil(SPI.skope.grad_free_dur/system.blockDurationRaster) * system.blockDurationRaster;
    if SPI.skope.grad_free_dur < SPI.skope.ttl_dur
        SPI.skope.grad_free_dur = SPI.skope.ttl_dur;
    end
    SPI.skope.delay    = mr.calcDuration(SPI.rf(1), SPI.gz_exc) + mr.calcDuration(SPI.gz_reph(1));
    SPI.skope.trig_out = mr.makeDigitalOutputPulse('ext1', 'duration', SPI.skope.ttl_dur, 'delay', SPI.skope.delay - SPI.skope.grad_free_dur - SPI.skope.ttl_dur);
    SPI.skope.delay    = mr.makeDelay(SPI.skope.delay);
else
    SPI.skope = [];
    SPI.skope.flag_onoff = false;
end

%% adapt timings for TRs and TEs

% create or import TR list
if ~isfield(SPI, 'TR')
    SPI.TR = zeros(SPI.NR,1);
end
if isfield(SPI, 'mrf_import')
    if isfield(SPI.mrf_import, 'TRs')
        if numel(SPI.mrf_import.TRs) == SPI.NR
            SPI.TR = SPI.mrf_import.TRs;
        else
            error('incorrect number of TRs');
        end
    end
else
    if numel(SPI.TR) < SPI.NR
        SPI.TR = SPI.TR * ones(SPI.NR,1);
    end
end
SPI.TR = SPI.TR(:);

% create or import TE list
if ~isfield(SPI, 'TE')
    SPI.TE = zeros(SPI.NR,1);
end
if isfield(SPI, 'mrf_import')
    if isfield(SPI.mrf_import, 'TEs')
        if numel(SPI.mrf_import.TEs) == SPI.NR
            SPI.TE = SPI.mrf_import.TEs;
        else
            error('incorrect number of TEs');
        end
    end
else
    if numel(SPI.TE) < SPI.NR
        SPI.TE = SPI.TE * ones(SPI.NR,1);
    end
end
SPI.TE = SPI.TE(:);

% calculate TE filling delays
temp_te_min  = mr.calcDuration(SPI.gz_exc)/2 + mr.calcDuration(SPI.gz_reph(1)) + system.blockDurationRaster;
temp_te_fill = SPI.TE - temp_te_min;
if sum(temp_te_fill<0)>0
    warning('negative TE filling delay! TE corrected!');
end
for j = 1:SPI.NR
    if temp_te_fill(j) >= system.blockDurationRaster
        SPI.TE_delay(j,1) = mr.makeDelay( round(temp_te_fill(j)/system.blockDurationRaster) * system.blockDurationRaster );
    else
        SPI.TE_delay(j,1) = mr.makeDelay( system.blockDurationRaster );
        SPI.TE(j)         = temp_te_min;
    end
end
clear temp_te_min temp_te_fill j;

% calculate TR filling delay
for j = 1:SPI.NR
    temp_tr_min(j,1) = mr.calcDuration(SPI.rf, SPI.gz_exc) + mr.calcDuration(SPI.gz_reph(1)) + mr.calcDuration(SPI.TE_delay(j)) + mr.calcDuration(SPI.gx, SPI.gy, SPI.gz, SPI.adc) + mr.calcDuration(SPI.gz_spoil) + system.gradRasterTime;
end
temp_tr_fill = SPI.TR - temp_tr_min;
if sum(temp_tr_fill<0)>0
    warning('negative TR filling delay! TR corrected!');
end
for j = 1:SPI.NR
    if temp_tr_fill(j) >= system.blockDurationRaster
        SPI.TR_delay(j,1) = mr.makeDelay( round(temp_tr_fill(j)/system.blockDurationRaster) * system.blockDurationRaster );
    else
        SPI.TR_delay(j,1) = mr.makeDelay( system.blockDurationRaster );
        SPI.TR(j)         = temp_tr_min(j);
    end
end
clear temp_tr_min temp_tr_fill j;

%% export reference k-space trajectory for reconstruction
if nargin<4
    flag_plot = 0;
end
if flag_plot==1 && strcmp(SPI.mode_2D_3D, '3D')
    flag_plot = 2;
end
temp_seq = mr.Sequence(system);
temp_seq.addBlock(SPI.rf(1));
temp_seq.addBlock(SPI.gx, SPI.gy, SPI.gz, SPI.adc, mr.makeRotation(0,0));
ktraj_ref = pulseq_get_ktraj(temp_seq, flag_plot);
clear temp_seq;

% note:
% this is how the k-space trajectory can be calculate for the other quaternions:
% ktraj_reco = repmat(ktraj_ref, 1, 1, SPI.proj.Nid);
% for j = 1:SPI.proj.Nid
%     ktraj_reco(:,:,j) = mr.aux.quat.toRotMat(SPI.rot(j).rotQuaternion) * ktraj_ref;
% end
% ktraj_reco = permute(ktraj_reco, [1,3,2]);

%% reorder struct field names
temp_names = fieldnames(SPI);
[~, idx]   = sort(lower(temp_names));
SPI        = orderfields(SPI, temp_names(idx));

end
