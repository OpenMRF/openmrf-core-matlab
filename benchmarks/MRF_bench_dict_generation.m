function out = MRF_bench_dict_generation(cfg)
%MRF_BENCH_DICT_GENERATION Benchmark the MRF dictionary forward simulation.
%
% out = MRF_BENCH_DICT_GENERATION(cfg) times the EPG (MRF_sim_EPG) and/or
% isochromat-based Bloch (MRF_sim_BLOCH) dictionary generation backends
% for the sequence in cfg.seq_file, across a sweep of dictionary sizes
% and in single-threaded and/or parallel-pool threading modes. Headline
% metric: per-atom generation time (ms/atom, atoms/s).
%
% cfg fields (all optional except seq_file):
%   seq_file      .seq file to simulate                        (required)
%   f0            [Hz] Larmor frequency                          (128e6)
%   dt            [s]  simulation raster time                     (1e-6)
%   backends      cell array, any of {'EPG','BLOCH'}            ({'EPG'})
%   N_iso         isochromats for the BLOCH backend                 (200)
%   N_dict_sweep  dictionary sizes for the scaling sweep    ([1 10 100 1000])
%   N_dict_sat    struct, large "saturated" dict size per backend name;
%                 timed in pool mode only, skipped for backends with no
%                 matching field                                    (-)
%   threading     'single' | 'pool' | 'both'                      ('both')
%   n_repeats     timed runs per point -> median                      (3)
%   warmup        discard one cold run first (parpool spin-up,      (true)
%                 cold pre-sim caches)
%   report_parse  also time & report .seq parse time                (true)
%   verbose       print progress/summary to the console              (true)
%
% out fields:
%   T         results table, one row per (backend, threading, N_dict)
%   R         same results as a struct array
%   env       host/run metadata (MATLAB version, core count, git rev, ...)
%   t_parse   median .seq parse time [s] (NaN if cfg.report_parse is false)
%   N_TR      number of TRs/ADCs in the parsed sequence
%   N_steps   number of raw operator steps in the parsed sequence
%   cfg       the cfg struct actually used (defaults filled in)
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany

cfg = fill_defaults(cfg);

assert(exist('MRF_read_seq_file','file')==2, ...
    'OpenMRF not on path -> run install_OpenMRF first.');

hasPCT = ~isempty(ver('parallel'));
if ~hasPCT
    if cfg.verbose
        warning('Parallel Computing Toolbox not found -> only single-threaded timing.');
    end
    cfg.threading = 'single';
end

env.datetime   = datetime('now');
env.matlab     = version;
env.computer   = computer;
env.numcores   = feature('numcores');
env.maxThreads = maxNumCompThreads;
try
    [s, g] = system('git rev-parse --short HEAD');
    if s==0, env.git = strtrim(g); else, env.git = ''; end
catch
    env.git = '';
end

if cfg.verbose
    fprintf('\n==================  MRF dictionary generation benchmark  ==================\n');
    fprintf('host: %d physical cores | MATLAB %s | %s\n', ...
            env.numcores, env.matlab, char(env.datetime));
end

%% ------------------------- parse sequence (once) ---------------------------

t_parse = NaN;
if cfg.report_parse
    t_parse = bench_median(@() MRF_read_seq_file(cfg.seq_file, cfg.f0, [],[],[],[],[],[], cfg.dt, 0), cfg.n_repeats);
end
[~, SIM0] = MRF_read_seq_file(cfg.seq_file, cfg.f0, [],[],[],[],[],[], cfg.dt, 0);
N_TR    = sum(SIM0.ID==0);
N_steps = numel(SIM0.ID);
if cfg.verbose
    fprintf('sequence: %s | TRs/ADCs: %d | raw operator steps: %d | parse time: %s\n', ...
            cfg.seq_file, N_TR, N_steps, fmt_t(t_parse));
end

%% ------------------- dictionary parameter pool + iso grid ------------------

N_sat_list = structfun(@(v) v, cfg.N_dict_sat);
N_max      = max([cfg.N_dict_sweep(:); N_sat_list(:)]);
P_pool     = make_param_pool(N_max);
z_iso      = linspace(-1/2, 1/2, cfg.N_iso)' * 2 * 8e-3;   % 8 mm slice, x2 oversampled
dw0_distr  = randn(cfg.N_iso, 1) * 10 * 2*pi;              % T2' isochromat distribution

%% ------------------------------- threading ---------------------------------

switch lower(cfg.threading)
    case 'single', th_modes = {'single'};
    case 'pool',   th_modes = {'pool'};
    case 'both',   th_modes = {'pool','single'};   % pool first: builds cold pre-sim caches in parallel
    otherwise,     error('cfg.threading must be single | pool | both');
end

%% --------------------------------- run -------------------------------------

R = struct('backend',{},'threading',{},'workers',{},'N_dict',{},'N_iso',{}, ...
           't_pre_s',{},'t_sim_s',{},'ms_per_atom',{},'atoms_per_s',{},'saturated',{});

for ib = 1:numel(cfg.backends)
    backend = cfg.backends{ib};

    for it = 1:numel(th_modes)
        th = th_modes{it};
        [restore_th, n_workers] = set_threading(th, env.numcores, hasPCT);

        % sweep points (+ one saturated point in pool mode, if configured)
        N_list = cfg.N_dict_sweep;
        sat    = false(size(N_list));
        if strcmp(th,'pool') && hasPCT && isfield(cfg.N_dict_sat, backend)
            N_list = [N_list, cfg.N_dict_sat.(backend)]; %#ok<AGROW>
            sat    = [sat, true];                        %#ok<AGROW>
        end

        % warm-up: full pre + sim at a small size -> discard (pays parpool spin-up
        % and the cold pre-simulation caches: adiabatic TM, slice profiles)
        if cfg.warmup
            Pw   = subset_P(P_pool, 8);
            SIMw = MRF_sim_pre(SIM0, Pw, z_iso, backend, 0, 0);
            quiet(@() run_backend(backend, SIMw, Pw, z_iso, dw0_distr));
            clear Pw SIMw;
        end

        if cfg.verbose
            fprintf('\n---- backend: %-5s | threading: %-6s | workers: %d ----\n', backend, th, n_workers);
            fprintf('%10s  %10s  %12s  %12s  %14s  %s\n', 'N_dict','t_pre[s]','t_sim[s]','ms/atom','atoms/s','');
        end
        for k = 1:numel(N_list)
            N    = N_list(k);
            Psub = subset_P(P_pool, N);

            t_pre = bench_median(@() MRF_sim_pre(SIM0, Psub, z_iso, backend, 0, 0), max(1,cfg.n_repeats-1));
            SIMp  = MRF_sim_pre(SIM0, Psub, z_iso, backend, 0, 0);
            t_sim = bench_median(@() quiet(@() run_backend(backend, SIMp, Psub, z_iso, dw0_distr)), cfg.n_repeats);

            ms_atom  = t_sim / N * 1e3;
            atoms_s  = N / t_sim;
            if cfg.verbose
                tag = ''; if sat(k), tag = '<- saturated'; end
                fprintf('%10d  %10.4f  %12.4f  %12.4f  %14.1f  %s\n', N, t_pre, t_sim, ms_atom, atoms_s, tag);
            end

            R(end+1) = struct('backend',backend,'threading',th,'workers',n_workers, ...
                              'N_dict',N,'N_iso',cfg.N_iso,'t_pre_s',t_pre,'t_sim_s',t_sim, ...
                              'ms_per_atom',ms_atom,'atoms_per_s',atoms_s,'saturated',sat(k)); %#ok<AGROW>
            clear SIMp Psub;
        end

        clear restore_th;   % onCleanup -> restores threading state
    end
end

T = struct2table(R);

%% ------------------------------- summary -----------------------------------

if cfg.verbose
    fprintf('\n========================  summary (per-atom)  ========================\n');
    for ib = 1:numel(cfg.backends)
        backend = cfg.backends{ib};
        for it = 1:numel(th_modes)
            th = th_modes{it};
            % representative point: largest non-saturated sweep entry
            sel = strcmp(T.backend,backend) & strcmp(T.threading,th) & ~T.saturated;
            if any(sel)
                sub = T(sel,:); [~,i] = max(sub.N_dict);
                fprintf('%-5s | %-6s | %5d-atom IR-FISP dict: %s  (%.3f ms/atom, %.0f atoms/s)\n', ...
                    backend, th, sub.N_dict(i), fmt_t(sub.t_sim_s(i)), sub.ms_per_atom(i), sub.atoms_per_s(i));
            end
            sel = strcmp(T.backend,backend) & strcmp(T.threading,th) & T.saturated;
            if any(sel)
                sub = T(sel,:);
                fprintf('%-5s | %-6s | saturated (%d atoms): %.4f ms/atom (%.0f atoms/s)\n', ...
                    backend, th, sub.N_dict(1), sub.ms_per_atom(1), sub.atoms_per_s(1));
            end
        end
    end
    fprintf('=====================================================================\n\n');
end

out.T       = T;
out.R       = R;
out.env     = env;
out.t_parse = t_parse;
out.N_TR    = N_TR;
out.N_steps = N_steps;
out.cfg     = cfg;

end

%% ============================ local functions =============================

function cfg = fill_defaults(cfg)
    assert(isfield(cfg,'seq_file') && ~isempty(cfg.seq_file), ...
        'cfg.seq_file is required.');
    if ~isfield(cfg,'f0'),            cfg.f0           = 128e6; end
    if ~isfield(cfg,'dt'),            cfg.dt           = 1e-6; end
    if ~isfield(cfg,'backends'),      cfg.backends     = {'EPG'}; end
    if ~isfield(cfg,'N_iso'),         cfg.N_iso        = 200; end
    if ~isfield(cfg,'N_dict_sweep'),  cfg.N_dict_sweep = [1 10 100 1000]; end
    if ~isfield(cfg,'N_dict_sat'),    cfg.N_dict_sat   = struct(); end
    if ~isfield(cfg,'threading'),     cfg.threading    = 'both'; end
    if ~isfield(cfg,'n_repeats'),     cfg.n_repeats    = 3; end
    if ~isfield(cfg,'warmup'),        cfg.warmup       = true; end
    if ~isfield(cfg,'report_parse'),  cfg.report_parse = true; end
    if ~isfield(cfg,'verbose'),       cfg.verbose      = true; end
end

function P = make_param_pool(N)
    % log-spaced T1/T2 grid clipped to T2<T1, subsampled (or tiled) to N atoms,
    % plus realistic B1+/B0 spreads so all simulation operators are exercised.
    q.T1.range = [0.01,  3.5]; q.T1.factor = 1.05;
    q.T2.range = [0.001, 1.5]; q.T2.factor = 1.05;
    q = MRF_get_param_dict(q, {'T2<T1'});
    ng  = numel(q.T1);
    if ng >= N
        idx = round(linspace(1, ng, N));
    else
        idx = mod(0:N-1, ng) + 1;          % tile if we need more atoms than grid points
    end
    P.T1  = q.T1(idx);
    P.T2  = q.T2(idx);
    P.db1 = 1 + 0.1 * linspace(-1, 1, N)';
    P.dw0 = 10 * 2*pi * linspace(-1, 1, N)';
end

function Ps = subset_P(P, n)
    % take n atoms spread evenly across the pool (duplicates allowed)
    np  = numel(P.T1);
    idx = round(linspace(1, np, n));
    fn  = fieldnames(P);
    for j = 1:numel(fn)
        v = P.(fn{j});
        Ps.(fn{j}) = v(idx);
    end
end

function M = run_backend(backend, SIMp, P, z, dw0)
    switch backend
        case 'EPG',   M = MRF_sim_EPG(SIMp, P, 0);
        case 'BLOCH', M = MRF_sim_BLOCH(SIMp, P, z, dw0, 0);
        otherwise,    error('unknown backend: %s', backend);
    end
end

function t = bench_median(fn, n)
    % median wall-clock time of n calls to fn (return values discarded)
    ts = zeros(n,1);
    for k = 1:n
        t0 = tic; fn(); ts(k) = toc(t0);
    end
    t = median(ts);
end

function quiet(fn)
    % run fn, suppressing its console output and discarding its (large) return value
    [~] = evalc('out_ = fn(); clear out_;');
end

function [restorer, n_workers] = set_threading(mode, ncores, hasPCT)
    % configure threading; returns an onCleanup object that restores prior state
    switch mode
        case 'single'
            prev.maxT = maxNumCompThreads(1);
            prev.pool = [];
            prev.auto = [];
            if hasPCT
                ps = parallel.Settings;
                prev.auto = ps.Pool.AutoCreate;
                ps.Pool.AutoCreate = false;
                prev.pool = gcp('nocreate');
                if ~isempty(prev.pool), delete(prev.pool); end
            end
            n_workers = 1;
            restorer  = onCleanup(@() restore_single(prev, hasPCT));

        case 'pool'
            % keep current thread count; ensure a pool exists
            maxNumCompThreads('automatic');
            if hasPCT
                ps = parallel.Settings; ps.Pool.AutoCreate = true;
                p  = gcp('nocreate');
                if isempty(p), p = parpool; end
                n_workers = p.NumWorkers;
            else
                n_workers = 1;
            end
            restorer = onCleanup(@() []);   % leave the pool running

        otherwise
            error('threading mode must be single or pool');
    end
end

function restore_single(prev, hasPCT)
    maxNumCompThreads(prev.maxT);
    if hasPCT && ~isempty(prev.auto)
        ps = parallel.Settings; ps.Pool.AutoCreate = prev.auto;
    end
end

function s = fmt_t(t)
    if isnan(t),      s = 'n/a';
    elseif t < 1,     s = sprintf('%.0f ms', t*1e3);
    elseif t < 90,    s = sprintf('%.2f s', t);
    else,             s = sprintf('%.1f min', t/60);
    end
end
