%% bench_mrf_generation.m
%
% Demo runner for MRF_bench_dict_generation: benchmarks the MRF dictionary
% forward simulation for a typical IR-FISP MRF sequence (1000-TR Yun-Jiang
% flip-angle pattern, bundled as test_mrf.seq), then saves the results
% table and plots ms/atom vs. dictionary size.
%
% Edit the cfg block below to change scope (drop BLOCH, shrink the sweep,
% etc). To benchmark a different sequence, backend set, or dictionary-size
% sweep programmatically, call MRF_bench_dict_generation(cfg) directly
% from your own script instead of editing this file.
%
% Usage:
%   bench_mrf_generation        % run with the defaults in the cfg block below
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.2, 07.07.2026

clear; clc;

%% ------------------------------ configuration ------------------------------

cfg.seq_file      = 'test_mrf.seq';          % bundled 1000-TR IR-FISP "yun" pattern
cfg.f0            = 128e6;                    % [Hz] Larmor frequency
cfg.dt            = 1e-6;                     % [s]  simulation raster time

cfg.backends      = {'EPG', 'BLOCH'};        % use {'EPG'} for the lightest run
cfg.N_iso         = 200;                      % isochromats for the BLOCH backend
cfg.N_dict_sweep  = [1 10 100 1000];    % dictionary sizes for the scaling sweep
cfg.N_dict_sat    = struct('EPG', 5e4, ...    % large "saturated" dict, per backend
                           'BLOCH', 1e4);     %   (pool mode only; BLOCH is slower)

cfg.threading     = 'both';                   % 'single' | 'pool' | 'both'
cfg.n_repeats     = 3;                        % timed runs per point -> median
cfg.warmup        = true;                     % discard one cold run first
                                              %   (parpool spin-up, cold pre-sim caches)

cfg.report_parse  = true;                     % also report .seq parse time
cfg.verbose       = true;                     % print progress/summary to the console

save_results      = true;                     % write <name>_<timestamp>.mat / .csv
make_plot         = true;                     % log-x plot: ms/atom vs N_dict

%% --------------------------------- run --------------------------------------

out = MRF_bench_dict_generation(cfg);

T = out.T; R = out.R; env = out.env; %#ok<NASGU>
t_parse = out.t_parse; N_TR = out.N_TR; N_steps = out.N_steps; %#ok<NASGU>

%% ------------------------------- save --------------------------------------

if save_results
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    base  = ['bench_mrf_generation_' stamp];
    save([base '.mat'], 'T', 'R', 'cfg', 'env', 't_parse', 'N_TR', 'N_steps');
    try, writetable(T, [base '.csv']); catch, end
    fprintf('saved: %s.mat / .csv\n', base);
end

%% ------------------------------- plot --------------------------------------

if make_plot
    figure('Name','MRF dictionary generation benchmark');
    hold on; styles = {'-o','--s',':^','-.d'};
    ci = 0; leg = {};
    th_modes = unique(T.threading, 'stable');
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
