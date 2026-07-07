% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 06.07.2026

%% add sequence objects for cardio or abdominal MRF

% init loop counters for contrast preparations
loop_SAT    = 1;
loop_INV    = 1;
loop_T2     = 1;
loop_SL     = 1;
loop_MLEV   = 1;
loop_ADIASL = 1;

% reset rf spoiling increment
loop_rf_inc = 0;  
   
for loop_MRF = 1 : MRF.n_segm

	% TRID segment label for GE scanners
    seq.addTRID(['CMRF_' num2str(loop_MRF) '_' MRF.enc_list{loop_MRF}]);
 
    % Trigger: R-Wave
    if ~isempty(TRIG_IN)
        seq.addBlock(TRIG_IN);   
    end

    % Soft Delay and Dynamic Delay: adjust cardiac phase or set individual recovery times
    if ~isempty(MRF.delay_soft)
        seq.addBlock(system.blockDurationRaster, MRF.delay_soft);
    end
    seq.addBlock(MRF.delay_dynamic(loop_MRF));   

    % MRF Preparations:
    MRF_add_preparation();

    % Spiral Readouts
    for loop_seg = 1:SPI.nr
        loop_NR = (loop_MRF-1)*SPI.nr + loop_seg;
        SPI_add();
    end

end