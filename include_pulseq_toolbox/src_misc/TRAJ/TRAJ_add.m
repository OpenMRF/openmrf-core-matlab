% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% add sequnece objects for trajectory scans

% ----- Duyn mode: -----
% Duyn JH, Yang Y, Frank JA, van der Veen JW.
% Simple correction method for k-space trajectory deviations in MRI.
% J Magn Reson. 1998 May;132(1):150-3.
% doi: 10.1006/jmre.1998.1396
% 1 -> positive gx gradient  in positive x-slice
% 2 -> reference phase       in positive x-slice
% 3 -> positive gy gradient  in positive y-slice
% 4 -> reference phase       in positive y-slice

% ----- Robison mode: -----
% Robison RK, Li Z, Wang D, Ooi MB, Pipe JG.
% Correction of B0 eddy current effects in spiral MRI.
% Magn Reson Med. 2019 Apr;81(4):2501-2513.
% doi: 10.1002/mrm.27583
% 1 -> positive gx gradient  in positive x-slice
% 2 -> negative gx gradient  in positive x-slice
% 3 -> positive gx gradient  in negative x-slice
% 4 -> negative gx gradient  in negative x-slice
% 5 -> positive gy gradient  in positive y-slice
% 6 -> negative gy gradient  in positive y-slice
% 7 -> positive gy gradient  in negative y-slice
% 8 -> negative gy gradient  in negative y-slice

%% loop cases for Duyn or Robison mode

if strcmp(TRAJ.method, 'duyn')
    temp_rf = TRAJ.rf;    
    switch loop_traj
        case 1
            temp_g_slice = TRAJ.g_exc_x;
            temp_g_reph  = TRAJ.g_reph_x;
            temp_g_read  = TRAJ.gx(loop_NR);
        case 2
            temp_g_slice = TRAJ.g_exc_x;
            temp_g_reph  = TRAJ.g_reph_x;
            temp_g_read  = TRAJ.ref_delay;
        case 3
            temp_g_slice = TRAJ.g_exc_y;
            temp_g_reph  = TRAJ.g_reph_y;
            temp_g_read  = TRAJ.gy(loop_NR);
        case 4
            temp_g_slice = TRAJ.g_exc_y;
            temp_g_reph  = TRAJ.g_reph_y;
            temp_g_read  = TRAJ.ref_delay;
    end
end

if strcmp(TRAJ.method, 'robison')
    temp_rf_pos = TRAJ.rf;
    temp_rf_neg = TRAJ.rf;
    temp_rf_neg.freqOffset = -temp_rf_neg.freqOffset;
    switch loop_traj
        case 1
            temp_rf      = temp_rf_pos;
            temp_g_slice = TRAJ.g_exc_x;
            temp_g_reph  = TRAJ.g_reph_x;
            temp_g_read  = TRAJ.gx(loop_NR);
        case 2
            temp_rf      = temp_rf_pos;
            temp_g_slice = TRAJ.g_exc_x;
            temp_g_reph  = TRAJ.g_reph_x;
            temp_g_read  = TRAJ.gx(loop_NR);
            temp_g_read.waveform = -temp_g_read.waveform;
        case 3
            temp_rf      = temp_rf_neg;
            temp_g_slice = TRAJ.g_exc_x;
            temp_g_reph  = TRAJ.g_reph_x;
            temp_g_read  = TRAJ.gx(loop_NR);
        case 4
            temp_rf      = temp_rf_neg;
            temp_g_slice = TRAJ.g_exc_x;
            temp_g_reph  = TRAJ.g_reph_x;
            temp_g_read  = TRAJ.gx(loop_NR);
            temp_g_read.waveform = -temp_g_read.waveform;
        case 5
            temp_rf      = temp_rf_pos;
            temp_g_slice = TRAJ.g_exc_y;
            temp_g_reph  = TRAJ.g_reph_y;
            temp_g_read  = TRAJ.gy(loop_NR);
        case 6
            temp_rf      = temp_rf_pos;
            temp_g_slice = TRAJ.g_exc_y;
            temp_g_reph  = TRAJ.g_reph_y;
            temp_g_read  = TRAJ.gy(loop_NR);
            temp_g_read.waveform = -temp_g_read.waveform;
        case 7
            temp_rf      = temp_rf_neg;
            temp_g_slice = TRAJ.g_exc_y;
            temp_g_reph  = TRAJ.g_reph_y;
            temp_g_read  = TRAJ.gy(loop_NR);
        case 8
            temp_rf      = temp_rf_neg;
            temp_g_slice = TRAJ.g_exc_y;
            temp_g_reph  = TRAJ.g_reph_y;
            temp_g_read  = TRAJ.gy(loop_NR);
            temp_g_read.waveform = -temp_g_read.waveform;            
    end
end

%% rf spoiling
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
temp_adc             = TRAJ.adc;
temp_phase           = 117*pi/180 * loop_rf_inc^2;
temp_rf.phaseOffset  = mod(temp_phase,        2*pi);
temp_adc.phaseOffset = mod(temp_phase + pi/2, 2*pi);
loop_rf_inc          = loop_rf_inc + 1;

%% prevent adcDeadTime violation
temp_n_zero = round((mr.calcDuration(temp_adc)-mr.calcDuration(temp_g_read))/system.gradRasterTime) + 2;
if temp_n_zero > 0
    for j = 1:temp_n_zero
        temp_g_read.waveform  = [temp_g_read.waveform; 0];
        temp_g_read.tt        = [temp_g_read.tt; temp_g_read.tt(end)+system.gradRasterTime];
        temp_g_read.shape_dur = temp_g_read.shape_dur + system.gradRasterTime;
    end
end
clear j temp_n_zero;

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
if strcmp(temp_g_slice.channel, 'x')
    seq.addTRID('exc_x');
elseif strcmp(temp_g_slice.channel, 'y')
    seq.addTRID('exc_y');
end
seq.addBlock(TRAJ.Trec);
seq.addBlock(temp_rf, temp_g_slice);
seq.addBlock(temp_g_reph);
seq.addBlock(TRAJ.d1);

if loop_av > 0
    if strcmp(temp_g_read.channel, 'x')
        seq.addTRID(['read_x_' num2str(loop_NR) num2str(mod(loop_traj, 2))]);
    elseif strcmp(temp_g_read.channel, 'y')
        seq.addTRID(['read_y_' num2str(loop_NR) num2str(mod(loop_traj, 2))]);
    end
    seq.addBlock(temp_g_read, temp_adc);
else
    if strcmp(temp_g_read.channel, 'x')
        seq.addTRID(['dummy_x_' num2str(loop_NR) num2str(mod(loop_traj, 2))]);
    elseif strcmp(temp_g_read.channel, 'y')
        seq.addTRID(['dummy_y_' num2str(loop_NR) num2str(mod(loop_traj, 2))]);
    end
    seq.addBlock(temp_g_read);
end

seq.addTRID('spoil_z');
seq.addBlock(TRAJ.d1);
seq.addBlock(TRAJ.g_spoil_z);

%% clear temp objects
clear temp_rf temp_rf_pos temp_rf_neg temp_g_slice temp_g_reph temp_g_read temp_adc temp_phase;