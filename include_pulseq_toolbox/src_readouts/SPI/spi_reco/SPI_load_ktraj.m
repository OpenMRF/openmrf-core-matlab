function ktraj_reco = SPI_load_ktraj(PULSEQ)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    
    % This function converts the single 3D reference k-space trajectory to the
    % individual trajectory for all unique readouts. The function
    % is necessary for the pulseq version v1.5.1 which uses rotation
    % objects and Quaternions instead of mr.rotate() functions.
    % There are different cases how to load ktraj_reco from the PULSEQ backup to ensure backwards compatibility. 

    % ----- input: -----
    % PULSEQ.SPI:        [struct]         contains all SPI params and rotation objects for unique projections
    % PULSEQ.ktraj_ref:  [3 x Nadc]       3D reference k-space trajectory with only one single projection
    % PULSEQ.ktraj_reco: [2 x Nid x Nadc] 2D reference k-space trajectory with all unique projections

    % ----- output: -----
    % ktraj_reco: [3 x Nid x Nadc] 3D k-space trajectory for all unique projections

    if isfield(PULSEQ.SPI, 'rot') % pulseq v.1.5.1
        ktraj_reco = zeros(size(PULSEQ.ktraj_ref,1), size(PULSEQ.ktraj_ref,2), numel(PULSEQ.SPI.rot));
        for j = 1:numel(PULSEQ.SPI.rot)
            ktraj_reco(:,:,j) = mr.aux.quat.toRotMat(PULSEQ.SPI.rot(j).rotQuaternion) * PULSEQ.ktraj_ref;
        end
        ktraj_reco = permute(ktraj_reco, [1,3,2]);
    else % pulseq v.1.5.0 or older
        ktraj_reco = PULSEQ.ktraj_reco;       
    end            

end