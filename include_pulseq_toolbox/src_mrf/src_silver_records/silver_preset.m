function spec = silver_preset(name)
% silver_preset
% Resolve a named "record-flavour" preset into a concrete build spec for
% make_silver_record. Presets are small JSON files in ./presets/ recalled by a
% memorable adj-adj-noun codename (e.g. 'calm-silver-heron') or a short alias
% (e.g. 'default'). The pattern is referenced by name and resolved via
% MRF_get_FAs_TRs; tissues via default_silver_tissues or MRF_get_param_dict.
%
% ---------------------------------------------------------------- inputs ---
% name : preset codename or alias (char/string, case-insensitive).
%
% --------------------------------------------------------------- outputs ---
% spec : struct ready for make_silver_record, with fields
%        .name .title         preset identity
%        .P .FAs .TRs .opt     resolved simulation inputs
%        .meta                 metadata for make_silver_record (stamps the
%                              record's .preset block: name/title/spec_sha256)
%        .spec_sha256          hash of the numeric excitation recipe
%                              (FAs, TRs, TE, T1, T2) -- language-neutral, so a
%                              sibling tool that reads the same values reproduces
%                              it byte-for-byte
%        .def_file             path of the resolved JSON definition
%
% Usage:
%   spec = silver_preset('default');
%   rec  = make_silver_record(spec.P, spec.FAs, spec.TRs, spec.opt, spec.meta);
% or in one step: make_silver_record_from_preset('default').
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    [def, file] = load_preset_def(name);

    % ----- flip-angle / TR pattern -----
    switch lower(def.pattern.source)
        case 'named'
            [FAs, TRs] = MRF_get_FAs_TRs(def.pattern.name, 0);
        otherwise
            error('silver_preset:pattern', 'unsupported pattern source: %s', def.pattern.source);
    end
    FAs = FAs(:);
    TRs = TRs(:);

    % ----- tissue table -----
    switch lower(def.tissues.mode)
        case 'named'
            names = {};
            if isfield(def.tissues, 'names') && ~isempty(def.tissues.names)
                names = cellstr(def.tissues.names);
            end
            [P, labels] = default_silver_tissues(names);
        case 'grid'
            g = def.tissues.grid;
            q.T1 = grid_axis(g.T1);
            q.T2 = grid_axis(g.T2);
            constraints = {};
            if isfield(g, 'constraints') && ~isempty(g.constraints)
                constraints = cellstr(g.constraints);
            end
            q = MRF_get_param_dict(q, constraints);
            P.T1 = q.T1(:); P.T2 = q.T2(:); labels = {};
        otherwise
            error('silver_preset:tissues', 'unknown tissue mode: %s', def.tissues.mode);
    end

    % ----- sequence options -----
    opt = struct('TE',0, 'spoil_twists',1, 'rf_phase',0, 'ref_phase',0, ...
                 't2_nrefocus',2, 'prep',{{}});
    if isfield(def, 'options') && ~isempty(def.options)
        o  = def.options;
        fn = fieldnames(o);
        for i = 1:numel(fn)
            f = fn{i};
            if strcmpi(f, 'prep')
                opt.prep = normalize_prep(o.prep);
            else
                opt.(f) = o.(f);
            end
        end
    end

    % ----- fingerprint of the numeric excitation recipe -----
    TEvec = opt.TE(:); if isscalar(TEvec), TEvec = TEvec*ones(numel(FAs),1); end
    v = double([real(FAs).', imag(FAs).', TRs.', TEvec.', P.T1(:).', P.T2(:).']);
    spec_sha256 = silver_sha256(typecast(v, 'uint8'));

    % ----- metadata that stamps the produced record -----
    meta.tissue_mode   = lower(def.tissues.mode);
    meta.tissue_labels = labels;
    meta.note          = def.title;
    meta.created_by    = sprintf('silver_preset:%s', def.name);
    meta.preset        = struct('name', def.name, 'title', def.title, ...
                                'alias', getf(def, 'alias', ''), ...
                                'spec_sha256', spec_sha256);

    spec.name        = def.name;
    spec.title       = def.title;
    spec.P           = P;
    spec.FAs         = FAs;
    spec.TRs         = TRs;
    spec.opt         = opt;
    spec.meta        = meta;
    spec.spec_sha256 = spec_sha256;
    spec.def_file    = file;
end

% ------------------------------------------------------------- helpers ---
function [def, file] = load_preset_def(name)
    d = fullfile(fileparts(mfilename('fullpath')), 'presets');
    files = dir(fullfile(d, '*.json'));
    if isempty(files)
        error('silver_preset:nopresets', 'no preset .json files found in %s', d);
    end
    name  = char(name);
    avail = {};
    for k = 1:numel(files)
        f   = fullfile(d, files(k).name);
        cur = jsondecode(fileread(f));
        nm  = getf(cur, 'name',  '');
        al  = getf(cur, 'alias', '');
        avail{end+1} = nm; %#ok<AGROW>
        if strcmpi(nm, name) || (~isempty(al) && strcmpi(al, name))
            def = cur; file = f; return;
        end
    end
    error('silver_preset:unknown', 'unknown preset "%s". Available: %s', ...
          name, strjoin(avail, ', '));
end

function ax = grid_axis(a)
    % build an MRF_get_param_dict axis spec from a JSON {range, factor|step}
    ax.range = a.range(:).';
    if isfield(a, 'factor') && ~isempty(a.factor)
        ax.factor = a.factor;
    elseif isfield(a, 'step') && ~isempty(a.step)
        ax.step = a.step;
    else
        error('silver_preset:grid', 'grid axis needs a factor or a step.');
    end
end

function prep = normalize_prep(raw)
    % Coerce jsondecode output into MRF_sim_irfisp_epg's prep format:
    % a row cell of blocks, each block a row cell {type(char), val(double)[, at]}.
    if isempty(raw), prep = {}; return; end
    if ~iscell(raw), raw = {raw}; end
    % a single block ["inv",0.02] may decode to a cell whose 1st entry is the
    % type string -> wrap it so it is treated as one block, not a list.
    if ~isempty(raw) && (ischar(raw{1}) || isstring(raw{1}))
        raw = {raw};
    end
    prep = cell(1, numel(raw));
    for i = 1:numel(raw)
        blk = raw{i};
        if ~iscell(blk), blk = num2cell(blk); end
        blk = blk(:).';
        blk{1} = char(blk{1});
        for j = 2:numel(blk), blk{j} = double(blk{j}); end
        prep{i} = blk;
    end
end

function v = getf(s, f, dflt)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
