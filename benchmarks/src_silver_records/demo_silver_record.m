%% demo_silver_record.m
%
% Driver: build a "silver-dictionary record" for the ideal IR-FISP (+optional
% T2-prep) EPG signal generator and write it out as a self-contained JSON file
% for cross-tool validation.
%
% A silver record bundles, in one interoperable file:
%   - the flip-angle pattern and TR pattern used,
%   - the sequence options (TE, spoiler, preparations),
%   - the (T1,T2) tissue table,
%   - the resulting complex signal dictionary Mxy,
%   - a layout manifest (shapes, dim order, units, complex & RF conventions),
%   - provenance (timestamp, git hash, machine id, MATLAB version),
%   - a payload checksum.
% A sibling simulator (e.g. a Python EPG port) exports the same structure and
% compare_silver_records() gauges whether the two agree.
%
% Requires OpenMRF on the path (run install_OpenMRF once).
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

clear; clc;

assert(exist('MRF_sim_irfisp_epg','file')==2, ...
    'MRF_sim_irfisp_epg not on path -> run install_OpenMRF first.');

%% -------- tissue set: 'named' (curated) or 'grid' (log-spaced) --------
tissue_mode = 'named';                      % 'named' | 'grid'
switch tissue_mode
    case 'named'
        [P, labels] = default_silver_tissues();
    case 'grid'
        q.T1.range = [0.20, 3.0]; q.T1.factor = 1.6;   % coarse: keeps N_dict small
        q.T2.range = [0.02, 0.6]; q.T2.factor = 1.6;
        q = MRF_get_param_dict(q, {'T2<T1'});
        P.T1 = q.T1; P.T2 = q.T2; labels = {};
    otherwise
        error('tissue_mode must be ''named'' or ''grid''.');
end
fprintf('tissue mode: %s  (%d atoms)\n', tissue_mode, numel(P.T1));

%% -------- flip-angle / TR pattern: canonical Yun-Jiang --------
[FAs, TRs] = MRF_get_FAs_TRs('yun', 0);      % [rad], [s]
NR = numel(FAs);

%% -------- sequence options (IR-FISP) --------
opt = struct();
opt.TE           = 2e-3;                      % [s] echo time
opt.spoil_twists = 1;                         % FISP unbalanced gradient spoiler
opt.t2_nrefocus  = 2;                         % (only used if a t2 prep is present)
opt.prep         = { {'inv', 20e-3} };        % inversion recovery at train start

%% -------- self-check: trust the generator before recording it --------
% T2-prep weighting on infinite T1 must equal exp(-tau/T2).
tau = 60e-3;
Pc.T1 = 1e6; Pc.T2 = 50e-3;
oc.TE = 0; oc.spoil_twists = 1; oc.t2_nrefocus = 2; oc.prep = {{'t2', tau}};
Mc  = MRF_sim_irfisp_epg(pi/2, 1, Pc, oc);
err = abs(abs(Mc(1)) - exp(-tau/Pc.T2));
fprintf('self-check T2-prep: err = %.2e  ->  %s\n', err, ternary(err<1e-3,'PASS','CHECK'));

%% -------- build the record --------
meta = struct();
meta.tissue_mode   = tissue_mode;
meta.tissue_labels = labels;
meta.note          = 'OpenMRF silver reference record (ideal IR-FISP EPG)';
meta.created_by    = 'demo_silver_record.m';

rec = make_silver_record(P, FAs, TRs, opt, meta);
fprintf('built record: %d readouts x %d tissues | payload sha256 %s...\n', ...
        rec.manifest.dims.shape(1), rec.manifest.dims.shape(2), ...
        rec.integrity.payload_sha256(1:12));

%% -------- write to <repo>/silver_records/ --------
here   = fileparts(mfilename('fullpath'));                 % .../src_mrf/src_silver_records
repo   = fileparts(fileparts(fileparts(here)));            % repo root
outdir = fullfile(repo, 'silver_records');
stamp  = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outfile = fullfile(outdir, sprintf('silver_irfisp_%s_%s.json', tissue_mode, stamp));
write_silver_record(rec, outfile);

%% -------- round-trip sanity: read back and self-compare --------
rec2 = read_silver_record(outfile);
compare_silver_records(rec, rec2, 1e-9);     % must be PASS / identical

%% -------- quick look at the recorded fingerprints --------
figure('Name', sprintf('silver record: %s', tissue_mode));
subplot(2,1,1);
plot(abs(FAs)*180/pi); grid on; ylabel('FA [deg]');
title(sprintf('Yun IR-FISP, NR=%d, %s tissues', NR, tissue_mode));
subplot(2,1,2);
plot(abs(rec2.Mxy)); grid on; xlabel('readout #'); ylabel('|Mxy|');
if ~isempty(labels), legend(labels, 'Interpreter','none', 'Location','northeast'); end

% ------------------------------------------------------------- helpers ---
function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
