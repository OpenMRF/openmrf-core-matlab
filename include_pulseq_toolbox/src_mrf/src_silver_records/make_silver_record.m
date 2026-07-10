function rec = make_silver_record(P, FAs, TRs, opt, meta)
% make_silver_record
% Build an in-memory "silver record": a self-describing bundle of a small
% IR-FISP (+ optional T2-prep) signal dictionary together with the flip-angle
% pattern, TR pattern, sequence options, tissue table, layout manifest and
% provenance. It is the interoperable reference a sibling simulator (e.g. a
% Python EPG port) is checked against for commensurability.
%
% This runs the ideal EPG generator MRF_sim_irfisp_epg and packages its result;
% it performs no file I/O. Use write_silver_record to serialize to JSON.
%
% ---------------------------------------------------------------- inputs ---
% P    : tissue parameter struct for MRF_sim_irfisp_epg. Requires .T1, .T2
%        (N_dict x 1) [s]. Optional: .dw0 [rad/s], .db1, .T1p, .m1p, ...
% FAs  : (NR x 1) flip angles [rad] (may be complex: |FAs|=flip, angle=phase).
% TRs  : (NR x 1) repetition times [s] (scalar -> constant TR).
% opt  : (optional) sequence options passed to MRF_sim_irfisp_epg
%        (.TE, .spoil_twists, .rf_phase, .prep, .ref_phase, .t2_nrefocus).
% meta : (optional) record metadata:
%        .tissue_mode   'named' | 'grid' | 'unspecified'   (default 'unspecified')
%        .tissue_labels {N_dict x 1} cellstr of tissue names (default {})
%        .note          free-text note stored in provenance
%        .created_by    free-text author/tool tag stored in provenance
%
% --------------------------------------------------------------- outputs ---
% rec : silver record struct with fields
%        .schema .sequence .manifest .inputs .outputs .provenance .integrity
%       ready to hand to write_silver_record.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    if nargin < 4 || isempty(opt),  opt  = struct(); end
    if nargin < 5 || isempty(meta), meta = struct(); end

    % ----- run the ideal EPG signal generator (the single source of truth) -----
    Mxy = MRF_sim_irfisp_epg(FAs, TRs, P, opt);

    % ----- resolve effective option defaults (mirror MRF_sim_irfisp_epg) -----
    def = struct('TE',0, 'spoil_twists',1, 'rf_phase',0, 'prep',{{}}, ...
                 'ref_phase',0, 't2_nrefocus',2);
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opt, fn{i}) || isempty(opt.(fn{i}))
            opt.(fn{i}) = def.(fn{i});
        end
    end

    % ----- normalise per-readout inputs for storage -----
    FAs = FAs(:);
    NR  = numel(FAs);
    N_dict = numel(P.T1);
    TRs = TRs(:); if isscalar(TRs), TRs = TRs*ones(NR,1); end
    TE  = opt.TE(:);  if isscalar(TE), TE = TE*ones(NR,1); end

    % ===================== schema / sequence =====================
    rec.schema   = 'openmrf.silver/1';
    rec.sequence = 'irfisp_t2prep_ideal';

    % preset stamp (when built from a named record-flavour preset)
    if isfield(meta, 'preset') && ~isempty(meta.preset)
        rec.preset = meta.preset;
    end

    % ===================== manifest (the commensurability contract) =====================
    man.dims.Mxy   = {'readout','tissue'};
    man.dims.shape = [NR, N_dict];
    man.dims.order = 'readout-major: rows = readouts, columns = tissue atoms';
    man.complex.stored_as  = 'real+imag';
    man.complex.convention = 'cartesian';
    man.rf.handedness = 'left-handed (alpha -> -alpha in MRF_sim_EPG/sim_rf)';
    man.rf.flip_unit  = 'rad';
    man.rf.notes      = '|RF| = flip angle, angle(RF) = rotation-axis phase';
    man.signal.M0             = 1.0;
    man.signal.demodulated_by = 'excitation phase (per readout), see MRF_sim_irfisp_epg';
    man.units = struct('FA','rad','TR','s','TE','s','T1','s','T2','s', ...
                       'dw0','rad/s','db1','dimensionless');
    rec.manifest = man;

    % ===================== inputs (fully regenerate the run) =====================
    in.FAs_real = real(FAs);
    in.FAs_imag = imag(FAs);
    in.TRs      = TRs;
    in.TE       = TE;
    in.spoil_twists = opt.spoil_twists;
    in.rf_phase     = opt.rf_phase;
    in.ref_phase    = opt.ref_phase;
    in.t2_nrefocus  = opt.t2_nrefocus;
    in.prep         = opt.prep;

    in.P.T1 = P.T1(:);
    in.P.T2 = P.T2(:);
    opt_pfields = {'dw0','db1','T1p','m1p','T1p_adia','ADC'};
    for k = 1:numel(opt_pfields)
        f = opt_pfields{k};
        if isfield(P, f) && ~isempty(P.(f))
            v = P.(f)(:);
            if isscalar(v), v = v*ones(N_dict,1); end
            in.P.(f) = v;
        end
    end

    in.tissue_mode   = getfield_default(meta, 'tissue_mode', 'unspecified');
    in.tissue_labels = getfield_default(meta, 'tissue_labels', {});
    rec.inputs = in;

    % ===================== outputs =====================
    out.Mxy_real = real(Mxy);
    out.Mxy_imag = imag(Mxy);
    rec.outputs  = out;

    % ===================== provenance =====================
    rec.provenance = build_provenance(meta);

    % ===================== integrity (numeric payload only) =====================
    rec.integrity.payload_sha256 = silver_checksum(rec);
end

% ------------------------------------------------------------- helpers ---
function v = getfield_default(s, f, dflt)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function pv = build_provenance(meta)
    pv.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    try
        pv.created_utc = char(datetime('now', 'TimeZone', 'UTC', ...
                                       'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
    catch
        pv.created_utc = '';
    end
    [pv.git_hash, pv.git_dirty] = git_state();
    host = get_hostname();
    pv.host           = host;
    pv.user           = get_username();
    pv.machine_id     = silver_sha256([host '|' computer]);
    pv.machine_id     = pv.machine_id(1:16);        % short pseudonymous id
    pv.os             = computer;
    pv.matlab_version = version;
    pv.cores          = feature('numcores');
    if isfield(meta, 'note'),       pv.note       = meta.note;       end
    if isfield(meta, 'created_by'), pv.created_by = meta.created_by; end
end

function [h, dirty] = git_state()
    h = ''; dirty = false;
    try
        [s, o] = system('git rev-parse --short HEAD');
        if s == 0, h = strtrim(o); end
        [s2, o2] = system('git status --porcelain');
        if s2 == 0, dirty = ~isempty(strtrim(o2)); end
    catch
    end
end

function h = get_hostname()
    h = getenv('COMPUTERNAME');
    if isempty(h), h = getenv('HOSTNAME'); end
    if isempty(h)
        try
            [s, o] = system('hostname');
            if s == 0, h = strtrim(o); end
        catch
        end
    end
    if isempty(h), h = 'unknown'; end
end

function u = get_username()
    u = getenv('USERNAME');
    if isempty(u), u = getenv('USER'); end
    if isempty(u), u = 'unknown'; end
end
