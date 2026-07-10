%% demo_silver_presets.m
%
% List the available record-flavour presets and build the default reference
% record from its codename. Presets are recalled by memorable adj-adj-noun
% names (or a short alias); the produced record is stamped with which preset
% and spec fingerprint made it.
%
% Requires OpenMRF on the path (run install_OpenMRF once).
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

clear; clc;

assert(exist('silver_preset','file')==2, ...
    'silver_preset not on path -> run install_OpenMRF first.');

%% -------- what is available --------
list_silver_presets();

%% -------- resolve one by alias (no simulation yet) --------
spec = silver_preset('default');
fprintf('resolved "%s" (%s)\n   NR=%d, N_dict=%d, spec %s...\n', ...
        spec.name, spec.title, numel(spec.FAs), numel(spec.P.T1), spec.spec_sha256(1:12));

%% -------- build + write the default reference record --------
[outfile, rec] = make_silver_record_from_preset('default');
fprintf('record preset stamp: %s / spec %s...\n', rec.preset.name, rec.preset.spec_sha256(1:12));

%% -------- other flavours are one call away --------
% [f2, rec2] = make_silver_record_from_preset('brisk-teal-marlin');   % IR + T2-prep
% [f3, rec3] = make_silver_record_from_preset('dense-amber-lattice'); % dense grid
