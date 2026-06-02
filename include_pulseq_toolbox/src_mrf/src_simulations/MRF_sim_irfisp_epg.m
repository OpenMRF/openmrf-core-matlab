function [Mxy, SIM] = MRF_sim_irfisp_epg(FAs, TRs, P, opt)
% MRF_sim_irfisp_epg
% Ideal (instantaneous-RF) EPG signal generator for a (T2-prepared) IR-FISP
% MRF sequence, built directly from {flip angles, TR pattern, preparations}.
%
% Purpose: isolate *signal generation* from all spatial/encoding machinery
% (no Pulseq, no FOV, no spiral, no .seq, no slice-profile / adiabatic
% pre-simulation). It hand-builds the operator list that MRF_sim_EPG walks
% and calls it directly. This is the cleanest 1:1 reference for validating a
% Python EPG implementation term-by-term.
%
% It is the lightweight analogue of the full pipeline
%   MRF_read_seq_file -> MRF_sim_pre -> MRF_sim_EPG
% collapsed to the ideal limit (instantaneous hard pulses, unit gradients).
%
% ---------------------------------------------------------------- inputs ---
% FAs : (NR x 1) flip angles [rad]. May be complex: |FAs| = flip angle,
%       angle(FAs) = RF excitation phase [rad] (e.g. for RF spoiling).
% TRs : (NR x 1) repetition times [s]  (scalar -> constant TR).
% P   : tissue parameter struct for MRF_sim_EPG. Requires .T1, .T2 (N_dict x 1).
%       Optional: .dw0 [rad/s] off-resonance, .db1 B1+ scaling, .T1p, ... .
% opt : (optional) struct of sequence options:
%   .TE           [s]  echo time, scalar or (NR x 1).            default 0
%   .spoil_twists [ ]  unit-gradient twists per TR (FISP).       default 1
%                       (must be >= 1; EPG needs >=1 dephasing state.
%                        N_states allocated by MRF_sim_EPG = ceil(sum|GZ|/2).)
%   .rf_phase     [rad] extra excitation phase (RF spoiling),
%                       scalar or (NR x 1).                       default 0
%   .prep         {cell} ordered preparation blocks (see below). By default
%                        they are placed BEFORE the readout train, but each
%                        block may carry a placement index to sit mid-train.
%                                                                  default {}
%   .ref_phase    [rad] reference phase for prep RF pulses.       default 0
%   .t2_nrefocus  [ ]  number of refocusing pulses in a T2 prep.  default 2
%
% Preparation blocks (opt.prep), an ordered cell array, each entry
% {type,val} or {type,val,at} where the optional 3rd element 'at' is the
% readout index the block is inserted immediately BEFORE:
%   {'inv', TI}      inversion (180 about x) + perfect crush + recovery TI [s]
%   {'t2',  tau}     T2 prep: 90 tipdown -> (MLEV 180s spanning tau) -> 90 tipup
%                    + perfect crush. tau is the T2-prep echo time [s].
%                    Ideal weighting on stored Mz is exp(-tau/T2).
%   {'delay', t}     pure free relaxation/recovery of t [s]
%   {'t2', tau, 500} same T2 prep, but emitted between readouts 499 and 500
%                    (the magnetization state carries through continuously).
% 'at' defaults to 1 (start of the train, the classic prep-before-readout
% layout); 'at' = NR+1 places the block after the final readout. To model
% segmented / interleaved schemes (e.g. cardiac MRF with a prep per heartbeat),
% give several blocks distinct 'at' indices, or splice operator lists.
%
% --------------------------------------------------------------- outputs ---
% Mxy : (NR x N_dict) complex transverse signal at each ADC (F0 state),
%       demodulated by the excitation phase.
% SIM : the operator-list struct passed to MRF_sim_EPG (for inspection /
%       export so the Python side can consume the identical program).
%
% Operator-list ID legend (see MRF_sim_EPG.m):
%   0  ADC (sample F0)     1  free relax + dephase (DT)     2  instantaneous RF
%   8  perfect crusher (zero transverse, keep Z)           10  gradient spoiler
%
% RF convention (MRF_sim_EPG > sim_rf): SIM.RF stores a complex flip angle,
%   |RF| = flip angle [rad], angle(RF) = rotation-axis phase [rad].
%   The internal rotation is left-handed (alpha -> -alpha in sim_rf). Mirror
%   exactly this convention in the Python reference. The demo self-check
%   validates the T2-prep weighting against exp(-tau/T2).
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 02.06.2026

    if nargin < 4 || isempty(opt), opt = struct(); end
    def = struct('TE',0, 'spoil_twists',1, 'rf_phase',0, 'prep',{{}}, ...
                 'ref_phase',0, 't2_nrefocus',2);
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opt, fn{i}) || isempty(opt.(fn{i}))
            opt.(fn{i}) = def.(fn{i});
        end
    end

    % ----- normalise per-readout inputs -----
    FAs = FAs(:);
    NR  = numel(FAs);
    TRs = TRs(:); if isscalar(TRs), TRs = TRs*ones(NR,1); end
    TE  = opt.TE(:);       if isscalar(TE),  TE  = TE *ones(NR,1); end
    ph  = opt.rf_phase(:); if isscalar(ph),  ph  = ph *ones(NR,1); end
    assert(numel(TRs)==NR && numel(TE)==NR && numel(ph)==NR, ...
        'TRs / TE / rf_phase must be scalar or (NR x 1).');
    assert(all(TRs >= TE - eps), 'each TR must be >= its TE.');
    assert(opt.spoil_twists >= 1, ...
        'spoil_twists must be >= 1 (EPG requires at least one dephasing state).');

    % ----- resolve prep placement (which readout each block sits before) -----
    nprep   = numel(opt.prep);
    prep_at = ones(nprep,1);                 % default: before the first readout
    for ip = 1:nprep
        if numel(opt.prep{ip}) >= 3 && ~isempty(opt.prep{ip}{3})
            prep_at(ip) = round(opt.prep{ip}{3});
        end
    end
    assert(all(prep_at >= 1 & prep_at <= NR+1), ...
        'prep placement index must be in [1, NR+1] (NR = %d).', NR);

    % event matrix columns: [ID, RF(complex), GZ, DT]
    E = zeros(0,4);

    % ===================== FISP readout train (preps interleaved) ==========
    PHI = zeros(NR,1);
    for n = 1:NR
        % emit any preparation block scheduled immediately before this readout
        for ip = find(prep_at(:).' == n)
            E = [E; build_prep(opt.prep{ip}, opt)];    %#ok<AGROW>
        end
        exc_phase = angle(FAs(n)) + ph(n);
        E = [E; rf_evt(abs(FAs(n)), exc_phase)];      %#ok<AGROW>  excitation
        E = [E; relax_evt(TE(n))];                     %#ok<AGROW>  -> TE
        E = [E; adc_evt()];                            %#ok<AGROW>  sample F0
        PHI(n) = exc_phase;                            %            demodulation
        E = [E; relax_evt(TRs(n)-TE(n))];              %#ok<AGROW>  -> end of TR
        E = [E; spoil_evt(opt.spoil_twists)];          %#ok<AGROW>  FISP spoiler
    end
    % preparation blocks scheduled after the final readout (at = NR+1)
    for ip = find(prep_at(:).' == NR+1)
        E = [E; build_prep(opt.prep{ip}, opt)];        %#ok<AGROW>
    end

    % ===================== assemble SIM struct =====================
    SIM.ID  = int32(real(E(:,1)));
    SIM.RF  = E(:,2);                 % complex flip angles
    SIM.GZ  = real(E(:,3));
    SIM.DB  = zeros(size(SIM.ID));    % no diffusion weighting
    SIM.DT  = real(E(:,4));
    SIM.PHI = PHI;
    SIM.TM  = [];                     % no adiabatic pulses in the ideal model
    SIM.seq_name = 'irfisp_t2prep_ideal';

    % ===================== simulate =====================
    Mxy = MRF_sim_EPG(SIM, P, 0);

end

% ------------------------------------------------------------- preparations ---
function E = build_prep(blk, opt)
% Emit the event matrix for one preparation block. blk is {type,val} or
% {type,val,at}; the optional 'at' placement index is consumed by the caller
% and ignored here.
    type = blk{1};
    val  = blk{2};
    switch lower(type)
        case {'inv','inversion'}
            E = [rf_evt(pi, opt.ref_phase);   ...   % 180 inversion
                 crush_evt();                 ...   % kill transverse
                 relax_evt(val)];                   % TI recovery
        case {'t2','t2prep'}
            E = t2prep_evts(val, opt.ref_phase, opt.t2_nrefocus);
        case {'delay','recovery','d'}
            E = relax_evt(val);
        otherwise
            error('MRF_sim_irfisp_epg:prep', 'unknown prep type: %s', type);
    end
end

% ----------------------------------------------------------------- events ---
function e = rf_evt(fa, phase),  e = [2, fa.*exp(1i*phase), 0, 0]; end
function e = adc_evt(),          e = [0, 0,                  0, 0]; end
function e = crush_evt(),        e = [8, 0,                  0, 0]; end
function e = spoil_evt(tw),      e = [10,0,                 tw, 0]; end
function e = relax_evt(dt)
    if dt > 0, e = [1, 0, 0, dt]; else, e = zeros(0,4); end
end

function E = t2prep_evts(tau, ref_phase, nref)
% Ideal T2 preparation: 90 tipdown -> nref MLEV 180s spread over tau -> 90 tipup
% -> perfect crush. Mirrors the block structure of T2_add.m in the ideal limit
% (instantaneous pulses, zero pulse duration). Net stored Mz weighting: exp(-tau/T2).
    if nargin < 3 || isempty(nref), nref = 2; end
    nref = max(1, round(nref));

    % MLEV-style alternating refocusing axes: +y, -y, +y, ...
    ref_axes = ref_phase + pi/2 * (-1).^(0:nref-1);

    % relaxation is split so the total transverse time equals tau:
    %   tau/(2*nref) at each end, tau/nref between successive refocusing pulses.
    t_end = tau/(2*nref);
    t_mid = tau/nref;

    E = rf_evt(pi/2, ref_phase + pi/2);          % 90 tipdown (about +y)
    E = [E; relax_evt(t_end)];
    for k = 1:nref
        E = [E; rf_evt(pi, ref_axes(k))];        % 180 refocusing
        if k < nref
            E = [E; relax_evt(t_mid)];
        end
    end
    E = [E; relax_evt(t_end)];
    E = [E; rf_evt(pi/2, ref_phase - pi/2)];     % 90 tipup (about -y = inverse tipdown)
    E = [E; crush_evt()];                        % crush residual transverse, keep Mz
end
