function [outpath, rec] = make_silver_record_from_preset(name, outdir)
% make_silver_record_from_preset
% One-step convenience: resolve a named preset, build its silver record and
% write it to disk.
%
% ---------------------------------------------------------------- inputs ---
% name   : preset codename or alias (see silver_preset / list_silver_presets).
% outdir : (optional) output directory.  default <repo>/silver_records
%
% --------------------------------------------------------------- outputs ---
% outpath : path of the written .json record.
% rec     : the in-memory record struct (carries rec.preset).
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    spec = silver_preset(name);
    rec  = make_silver_record(spec.P, spec.FAs, spec.TRs, spec.opt, spec.meta);

    if nargin < 2 || isempty(outdir)
        here   = fileparts(mfilename('fullpath'));            % .../src_silver_records
        repo   = fileparts(fileparts(fileparts(here)));       % repo root
        outdir = fullfile(repo, 'silver_records');
    end

    stamp   = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outpath = fullfile(outdir, sprintf('silver_%s_%s.json', spec.name, stamp));
    write_silver_record(rec, outpath);
end
