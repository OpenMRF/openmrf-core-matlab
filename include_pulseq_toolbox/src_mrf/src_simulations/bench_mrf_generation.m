%% bench_mrf_generation.m
%
% Light benchmark harness for the MRF dictionary forward simulation.
%
% Headline metric: per-atom dictionary generation time (ms/atom, atoms/s) for
% a typical IR-FISP MRF sequence (1000-TR Yun-Jiang flip-angle pattern, bundled
% as test_mrf.seq). Two backends are benchmarked:
%   'EPG'   -> MRF_sim_EPG   (extended phase graphs, instantaneous RF)
%   'BLOCH' -> MRF_sim_BLOCH (isochromat-based fast Bloch, slice-profile corrected)
% The slow brute-force Bloch simulation (MRF_sim_brute_force) is intentionally
% NOT benchmarked.
%
% Each backend is timed in two threading modes so results compare across systems:
%   'single' -> maxNumCompThreads(1), no parallel pool (parfor runs serially)
%   'pool'   -> default parallel pool (NumWorkers recorded with each result)
% In addition, a large "saturated" dictionary is timed in pool mode so that the
% per-atom estimate is no longer dominated by fixed parfor overhead.
%
% Peripheral datapoint: .seq file parsing time (MRF_read_seq_file).
%
% Usage:
%   bench_mrf_generation        % run with the defaults in the cfg block below
% Edit the cfg block to change scope (drop BLOCH, shrink the sweep, etc).
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 02.06.2026

clear; clc;

%% ------------------------------ configuration ------------------------------

cfg.seq_file      = 'test_mrf.seq';          % bundled 1000-TR IR-FISP "yun" pattern
cfg.f0            = 128e6;                    % [Hz] Larmor frequency
cfg.dt            = 1e-6;                     % [s]  simulation raster time

cfg.backends      = {'EPG', 'BLOCH'};        % use {'EPG'} for the lightest run
cfg.N_iso         = 200;                      % isochromats for the BLOCH backend
cfg.N_dict_sweep  = [1 10 100 1000 10000];    % dictionary sizes for the scaling sweep
cfg.N_dict_sat    = struct('EPG', 5e4, ...    % large "saturated" dict, per backend
                           'BLOCH', 1e4);     %   (pool mode only; BLOCH is slower)

cfg.threading     = 'both';                   % 'single' | 'pool' | 'both'
cfg.n_repeats     = 3;                        % timed runs per point -> median
cfg.warmup        = true;                     % discard one cold run first
                                              %   (parpool spin-up, cold pre-sim caches)

cfg.report_parse  = true;                     % also report .seq parse time
cfg.save_results  = true;                     % write <name>_<timestamp>.mat / .csv
cfg.make_plot     = true;                     % log-x plot: ms/atom vs N_dict

%% ------------------------------ environment --------------------------------

assert(exist('MRF_read_seq_file','file')==2, ...
    'OpenMRF not on path -> run install_OpenMRF first.');

hasPCT = ~isempty(ver('parallel'));
if ~hasPCT
    warning('Parallel Computing Toolbox not found -> only single-threaded timing.');
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

fprintf('\n==================  MRF dictionary generation benchmark  ==================\n');
fprintf('host: %d physical cores | MATLAB %s | %s\n', ...
        env.numcores, env.matlab, char(env.datetime));

%% ------------------------- parse sequence (once) ---------------------------

t_parse = NaN;
if cfg.report_parse
    t_parse = bench_median(@() MRF_read_seq_file(cfg.seq_file, cfg.f0, [],[],[],[],[],[], cfg.dt, 0), cfg.n_repeats);
end
[~, SIM0] = MRF_read_seq_file(cfg.seq_file, cfg.f0, [],[],[],[],[],[], cfg.dt, 0);
N_TR    = sum(SIM0.ID==0);
N_steps = numel(SIM0.ID);
fprintf('sequence: %s | TRs/ADCs: %d | raw operator steps: %d | parse time: %s\n', ...
        cfg.seq_file, N_TR, N_steps, fmt_t(t_parse));

%% ------------------- dictionary parameter pool + iso grid ------------------

N_max     = max([cfg.N_dict_sweep, cfg.N_dict_sat.EPG, cfg.N_dict_sat.BLOCH]);
P_pool    = make_param_pool(N_max);
z_iso     = linspace(-1/2, 1/2, cfg.N_iso)' * 2 * 8e-3;   % 8 mm slice, x2 oversampled
dw0_distr = randn(cfg.N_iso, 1) * 10 * 2*pi;              % T2' isochromat distribution

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

        % sweep points (+ one saturated point in pool mode)
        N_list = cfg.N_dict_sweep;
        sat    = false(size(N_list));
        if strcmp(th,'pool') && hasPCT
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

        fprintf('\n---- backend: %-5s | threading: %-6s | workers: %d ----\n', backend, th, n_workers);
        fprintf('%10s  %10s  %12s  %12s  %14s  %s\n', 'N_dict','t_pre[s]','t_sim[s]','ms/atom','atoms/s','');
        for k = 1:numel(N_list)
            N    = N_list(k);
            Psub = subset_P(P_pool, N);

            t_pre = bench_median(@() MRF_sim_pre(SIM0, Psub, z_iso, backend, 0, 0), max(1,cfg.n_repeats-1));
            SIMp  = MRF_sim_pre(SIM0, Psub, z_iso, backend, 0, 0);
            t_sim = bench_median(@() quiet(@() run_backend(backend, SIMp, Psub, z_iso, dw0_distr)), cfg.n_repeats);

            ms_atom  = t_sim / N * 1e3;
            atoms_s  = N / t_sim;
            tag      = ''; if sat(k), tag = '<- saturated'; end
            fprintf('%10d  %10.4f  %12.4f  %12.4f  %14.1f  %s\n', N, t_pre, t_sim, ms_atom, atoms_s, tag);

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

%% ------------------------------- save --------------------------------------

if cfg.save_results
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    base  = ['bench_mrf_generation_' stamp];
    save([base '.mat'], 'T', 'R', 'cfg', 'env', 't_parse', 'N_TR', 'N_steps');
    try, writetable(T, [base '.csv']); catch, end
    fprintf('saved: %s.mat / .csv\n', base);
end

%% ------------------------------- plot --------------------------------------

if cfg.make_plot
    figure('Name','MRF dictionary generation benchmark');
    hold on; styles = {'-o','--s',':^','-.d'};
    ci = 0; leg = {};
    for ib = 1:numel(cfg.backends)
        backend = cfg.backends{ib};
        for it = 1:numel(th_modes)
            th = th_modes{it}; ci = ci + 1;
            sel = strcmp(T.backend,backend) & strcmp(T.threading,th);
            sub = sortrows(T(sel,:),'N_dict');
            plot(sub.N_dict, sub.ms_per_atom, styles{mod(ci-1,numel(styles))+1}, 'LineWidth', 2, 'MarkerSize', 7);
            leg{end+1} = sprintf('%s / %s', backend, th); %#ok<AGROW>
        end
    end
    set(gca,'XScale','log','YScale','log','FontName','Arial','FontWeight','bold','LineWidth',2);
    grid on; xlabel('dictionary size N_{dict}'); ylabel('time per atom [ms]');
    title(sprintf('MRF dictionary generation: %s (%d TR)', cfg.seq_file, N_TR), 'Interpreter','none');
    legend(leg,'Location','northeast');
end

%% ============================ local functions =============================

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
