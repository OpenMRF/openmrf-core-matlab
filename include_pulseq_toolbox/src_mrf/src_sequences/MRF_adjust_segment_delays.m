% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 06.07.2026

%% add dummy MRF preparations

% init loop counters for contrast preparations
loop_SAT    = 1;
loop_INV    = 1;
loop_T2     = 1;
loop_SL     = 1;
loop_MLEV   = 1;
loop_ADIASL = 1;
  
for loop_MRF = 1 : MRF.n_segm

    % reset seq struct
    clear seq;
    seq = mr.Sequence(system);
 
    % MRF Preparations:
    MRF_add_preparation();

    % calc prep durations
    MRF.prep_durations(loop_MRF,1) = seq.duration;

end
clear seq;
seq = mr.Sequence(system);

%% adjust segment timings and trigger options

% MRF.prep_durations:     durations of the individual preparation modules
% MRF.acq_durations:      durations of the individual readout blocks
% MRF.prep_acq_durations: durations of the individual segments -> MRF.prep_durations + MRF.acq_durations
% MRF.delay_soft:         soft delay for shifting the acquisition window using the scanner GUI
% MRF.delay_dynamic:      individual recovery delays OR filling delay to balance the timings

% calc segment specific durations
MRF.acq_durations      = sum(reshape(SPI.TR, [MRF.nr, MRF.n_segm]))';
MRF.prep_acq_durations = MRF.prep_durations + MRF.acq_durations;

if ~isfield(MRF, 'mode_seg')
    MRF.mode_trig = '';
end

switch MRF.mode_seg
    case 'trigger' % this will enable the cardiac trigger; segment timings will be synchronized; soft delays will be used for adjusting the acq window on the scanner
        MRF.mode_trig    = 'on';
        MRF.seg_duration = [];
        MRF.rec_times    = [];
        TRIG_IN          = mr.makeTrigger('physio1', 'system', system, 'duration', system.blockDurationRaster); % ECG input trigger
        MRF.delay_soft   = mr.makeSoftDelay(0, 'acq_end', 'offset', -round(max(MRF.prep_acq_durations),4), 'factor', 1); % block_duration [s] = offset [s] + input [s] / factor
        if numel(unique(round(MRF.acq_durations,2))) > 1
            warning('different MRF.acq_durations detected! -> variable acq windows! -> check TRs!');
        end
        for loop_MRF = 1 : MRF.n_segm
            MRF.delay_dynamic(loop_MRF,1) = mr.makeDelay( round(( ...
                                            max(MRF.prep_acq_durations) - MRF.prep_acq_durations(loop_MRF) ...
                                            ) / system.gradRasterTime)*system.gradRasterTime + system.gradRasterTime );
        end
        MRF.min_RR = max(MRF.prep_acq_durations);
        disp(['minimum RR: ' num2str(round(MRF.min_RR*1e3)) 'ms']);

    case 'constant'
        MRF.mode_trig  = 'off';
        MRF.rec_times  = [];
        TRIG_IN        = [];
        MRF.delay_soft = [];
        for loop_MRF = 1 : MRF.n_segm
            MRF.delay_dynamic(loop_MRF,1) = mr.makeDelay( round(( ...
                                            MRF.seg_duration - MRF.prep_acq_durations(loop_MRF) ...
                                            ) / system.gradRasterTime)*system.gradRasterTime + system.gradRasterTime );
        end

    case 'variable'
        MRF.mode_trig    = 'off';
        MRF.seg_duration = [];
        TRIG_IN          = [];
        MRF.delay_soft   = [];
        for loop_MRF = 1 : MRF.n_segm
            MRF.delay_dynamic(loop_MRF,1) = mr.makeDelay( round(( ...
                                            MRF.rec_times(loop_MRF) ...
                                            ) / system.gradRasterTime)*system.gradRasterTime );
        end

    otherwise
        error('choose MRF.mode_seg: >trigger< or >constant< or >variable<');
end

%% delte loop counters
clear loop_SAT loop_INV loop_T2 loop_SL loop_MLEV loop_ADIASL loop_MRF;
