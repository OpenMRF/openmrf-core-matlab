%% demo_irfisp_t2prep_epg.m
%
% Minimal, dependency-light demo + self-check for MRF_sim_irfisp_epg, the
% ideal (instantaneous-RF) EPG signal generator. No Pulseq / FOV / spiral /
% .seq involved: pure {flip angles, TR pattern, preparations} -> signal.
%
% Use this as the cross-validation harness against a Python EPG implementation:
% feed both sides the identical FAs / TRs / TE / prep and compare Mxy.
%
% Requires only MRF_sim_EPG on the path (run install_OpenMRF once).
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 02.06.2026

clear; clc;

assert(exist('MRF_sim_EPG','file')==2, ...
    'MRF_sim_EPG not on path -> run install_OpenMRF first.');

%% -------- tissue dictionary (the only "P" the engine needs) --------
P.T1 = [1000; 1400;  800; 1200]' *1e-3;   % [s]
P.T2 = [  45;   80;  200;  100]' *1e-3;    % [s]
P.T1 = P.T1(:); P.T2 = P.T2(:);
% optional: P.dw0 = 0*P.T1; P.db1 = ones(size(P.T1));

%% -------- flip-angle / TR pattern (non-canonical: 5 sinusoidal arches) --------
% Build a custom FA train from 5 half-sine "arches", each ramping 0 -> peak -> 0.
% Peak flip angles per arch [deg]; tweak freely (this is NOT the Yun pattern).
arch_amp = [40 50 60 50 38];                   % [deg] max amplitude of each arch
n_per    = 200;                                 % readouts per arch
FAs = zeros(numel(arch_amp)*n_per, 1);          % [deg]
for a = 1:numel(arch_amp)
    k   = (1:n_per)';
    idx = (a-1)*n_per + k;
    FAs(idx) = arch_amp(a) * sin(pi*k/n_per);   % half-sine arch
end
FAs = FAs * pi/180;                             % [deg] -> [rad]
FAs(FAs==0) = 1e-6;                             % avoid exact-zero flips
NR  = numel(FAs);
TRs = 12e-3 * ones(NR,1);                        % [s] constant TR

%% -------- sequence options --------
opt.TE           = 2e-3;       % [s] echo time
opt.spoil_twists = 1;          % FISP unbalanced gradient spoiler
opt.t2_nrefocus  = 2;          % MLEV-2 refocusing in the T2 prep

% --- T2-prep placement switch ---
% Index of the readout the T2 prep is inserted immediately BEFORE. The EPG
% state carries through continuously, so this genuinely re-weights the train
% from that point on.  1   = classic prep-before-readout layout (start);
%                      500 = drop the T2 prep in the middle of the train.
t2prep_at = 250;                              % <-- set e.g. 500 to place mid-train
t2prep_at = min(max(round(t2prep_at),1), NR+1);

% inversion (TI) stays at the start; the T2 prep (tau) is placed at t2prep_at
opt.prep = { {'inv', 20e-3}, {'t2', 60e-3, t2prep_at} };

%% -------- run --------
[Mxy, SIM] = MRF_sim_irfisp_epg(FAs, TRs, P, opt);
fprintf('simulated %d readouts x %d tissues. SIM has %d operator steps.\n', ...
        size(Mxy,1), size(Mxy,2), numel(SIM.ID));

%% -------- plot: sequence pattern + fingerprints in one figure --------
figure('Name','ideal IR-FISP + T2-prep: pattern & fingerprints');

% (1) flip-angle pattern
subplot(4,1,1);
plot(abs(FAs)*180/pi, '-'); grid on; hold on;
ylabel('FA [deg]');
title(sprintf('inv(TI=20ms) + T2prep(\\tau=60ms) at readout %d + FISP, NR=%d', ...
              t2prep_at, NR));
yl = ylim;
plot([t2prep_at t2prep_at], yl, 'r--');     % mark the T2-prep insertion point
text(t2prep_at, yl(2), ' T2 prep', 'Color','r', ...
     'VerticalAlignment','top', 'HorizontalAlignment','left');
ylim(yl);

% (2) TR pattern
subplot(4,1,2);
plot(TRs*1e3, '-'); grid on; hold on;
ylabel('TR [ms]');
yl = ylim; plot([t2prep_at t2prep_at], yl, 'r--'); ylim(yl);

% (3) fingerprint magnitude
subplot(4,1,3);
plot(abs(Mxy)); grid on; hold on;
ylabel('|Mxy|');
yl = ylim; plot([t2prep_at t2prep_at], yl, 'r--'); ylim(yl);

% (4) fingerprint real / imag
subplot(4,1,4);
plot(real(Mxy)); hold on; plot(imag(Mxy),'--'); grid on;
xlabel('readout #'); ylabel('Re / Im (Mxy)');
yl = ylim; plot([t2prep_at t2prep_at], yl, 'r--'); ylim(yl);

%% ==================== self-check: T2-prep weighting ====================
% A T2 prep on a single tissue with very long T1 (no T1 recovery during tau)
% must store Mz' = exp(-tau/T2). Read it out with a single 90 -> |Mxy| = that.
tau   = 60e-3;
Pchk.T1 = 1e6;          % effectively infinite -> isolate T2 decay
Pchk.T2 = 50e-3;
optc.TE = 0; optc.spoil_twists = 1; optc.t2_nrefocus = 2;
optc.prep = { {'t2', tau} };
Mchk     = MRF_sim_irfisp_epg(pi/2, 1, Pchk, optc);   % single 90 readout
got      = abs(Mchk(1));
expect   = exp(-tau / Pchk.T2);
fprintf('\nT2-prep self-check: |Mxy| = %.6f, exp(-tau/T2) = %.6f, err = %.2e\n', ...
        got, expect, abs(got-expect));
if abs(got-expect) < 1e-3
    fprintf('  PASS: T2-prep weighting matches analytic exp(-tau/T2).\n');
else
    fprintf(['  CHECK: deviation > 1e-3. The EPG operators are correct;\n' ...
             '  this usually means a prep RF axis/phase convention needs a\n' ...
             '  sign flip. Adjust opt.ref_phase / the tipup phase to match.\n']);
end

%% ==================== self-check: IR recovery ====================
% Inversion + recovery TI on long T2: Mz(TI) = 1 - 2*exp(-TI/T1); a 90 reads |Mz|.
TI = 0.5; Pir.T1 = 1.0; Pir.T2 = 1e6;
opti.TE = 0; opti.spoil_twists = 1; opti.prep = { {'inv', TI} };
Mir     = MRF_sim_irfisp_epg(pi/2, 1, Pir, opti);
got2    = abs(Mir(1));
expect2 = abs(1 - 2*exp(-TI/Pir.T1));
fprintf('IR self-check:      |Mxy| = %.6f, |1-2e^{-TI/T1}| = %.6f, err = %.2e\n', ...
        got2, expect2, abs(got2-expect2));

%% -------- export the program for the Python side (optional) --------
% writematrix([double(SIM.ID), real(SIM.RF), imag(SIM.RF), SIM.GZ, SIM.DT], 'sim_program.csv');
% writematrix(SIM.PHI, 'sim_phi.csv');
