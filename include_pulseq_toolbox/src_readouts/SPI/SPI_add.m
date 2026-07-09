% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: spiral (SPI) -----------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- sequence loop counter: -----
% loop_rf_inc -> rf spoiling
% loop_kz     -> z partitions
% loop_NR     -> repetitions
% loop_skope  -> skope segments

%% init loop counters
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
if ~exist('loop_kz', 'var')
    loop_kz = 1;
end
if loop_NR<1
    temp_loop = 1;
else
    temp_loop = loop_NR;
end
if ~exist('loop_skope', 'var')
    loop_skope = [];
end

%% rf spoilng
if strcmp(SPI.spoil_rf_mode, 'import')
    temp_phase  = SPI.mrf_import.rf_phases(temp_loop);
else
    temp_phase  = SPI.spoil_rf_inc * loop_rf_inc ^ SPI.spoil_rf_pow;
    loop_rf_inc = loop_rf_inc + 1;
end
SPI.rf(temp_loop).phaseOffset = mod(temp_phase,        2*pi);
SPI.adc.phaseOffset           = mod(temp_phase + pi/2, 2*pi); % shift with pi/2: Mx -> real axis
clear temp_phase;

%% LIN label extension for United Imaging scanners
if flag_UI==1
    if loop_NR>0
        if ~exist('loop_lin', 'var')
            loop_lin = 1;
        end
        seq.addBlock(mr.makeLabel('SET','LIN', loop_lin));
        loop_lin = loop_lin + 1;
    end
end

%% add sequence blocks
if ~SPI.skope.flag_onoff
    seq.addTRID('spiral_slice_excitation');
    seq.addBlock(SPI.rf(temp_loop), SPI.gz_exc);
    seq.addBlock(SPI.gz_reph(loop_kz));
else
    if ismember(loop_NR, SPI.skope.TR_list(:, loop_skope))
        seq.addTRID('skope_acquisition');
        seq.addBlock(SPI.skope.delay, SPI.skope.trig_out);
    else
        seq.addTRID('skope_dummy');
        seq.addBlock(SPI.gz_exc);
        seq.addBlock(SPI.gz_reph(loop_kz));
    end
end
seq.addBlock(SPI.TE_delay(temp_loop));
if loop_NR<1
    seq.addTRID('spiral_dummy');
    seq.addBlock(SPI.gx, SPI.gy, SPI.gz);
else
    seq.addTRID('spiral_readout');
    if flag_rot==1
        seq.addBlock(SPI.rot(SPI.proj.id(temp_loop)), SPI.gx, SPI.gy, SPI.gz, SPI.adc); % requires interpreter v1.5.1
    else
        seq.addBlock(mr.rotate3D(SPI.rot(SPI.proj.id(temp_loop)).rotQuaternion, SPI.gx, SPI.gy, SPI.gz, SPI.adc)); % v1.5.0 and older
    end
end
seq.addBlock(SPI.gz_spoil(loop_kz));
seq.addBlock(SPI.TR_delay(temp_loop));
clear temp_loop