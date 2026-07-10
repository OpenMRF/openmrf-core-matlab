function rep = compare_silver_records(recA, recB, tol)
% compare_silver_records
% Gauge the commensurability of two silver records: first check that they are
% structurally comparable (schema / sequence / layout / units / conventions),
% then align tissue atoms by their (T1,T2) key and quantify the signal
% agreement (RMSE, max-abs deviation, complex correlation).
%
% This is the acceptance harness for a sibling simulator: export its dictionary
% as a silver record and compare against the MATLAB reference.
%
% ---------------------------------------------------------------- inputs ---
% recA, recB : silver record structs (from make_/read_silver_record) OR file
%              paths to .json records (loaded automatically).
% tol        : (optional) max-abs deviation treated as "matching".  default 1e-4
%
% --------------------------------------------------------------- outputs ---
% rep : report struct with fields
%        .compatible      logical, all structural checks passed
%        .messages        {cellstr} human-readable incompatibility notes
%        .n_A .n_B        atom counts
%        .n_matched       atoms in A found (by T1,T2) in B
%        .rmse .max_abs   global signal deviation over matched atoms
%        .min_correlation worst per-atom complex correlation
%        .within_tol      logical, max_abs <= tol
%        .per_atom        struct of per-atom rmse / max_abs / correlation
%      A summary is also printed to the console.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    if nargin < 3 || isempty(tol), tol = 1e-4; end
    if ischar(recA) || isstring(recA), recA = read_silver_record(recA); end
    if ischar(recB) || isstring(recB), recB = read_silver_record(recB); end

    msgs = {};
    ok   = true;

    % ----- 1) schema / sequence -----
    if ~strcmp(get_field(recA,'schema',''), get_field(recB,'schema',''))
        msgs{end+1} = sprintf('schema differs: "%s" vs "%s"', ...
            get_field(recA,'schema',''), get_field(recB,'schema','')); ok = false; %#ok<AGROW>
    end
    if ~strcmp(get_field(recA,'sequence',''), get_field(recB,'sequence',''))
        msgs{end+1} = sprintf('sequence differs: "%s" vs "%s"', ...
            get_field(recA,'sequence',''), get_field(recB,'sequence','')); ok = false; %#ok<AGROW>
    end

    % ----- 2) layout / conventions -----
    shA = double(recA.manifest.dims.shape(:).');
    shB = double(recB.manifest.dims.shape(:).');
    if ~isequal(shA(1), shB(1))
        msgs{end+1} = sprintf('readout count differs: %d vs %d', shA(1), shB(1)); ok = false; %#ok<AGROW>
    end
    if ~isequal(recA.manifest.complex, recB.manifest.complex)
        msgs{end+1} = 'complex storage convention differs'; ok = false; %#ok<AGROW>
    end
    if ~strcmp(recA.manifest.rf.handedness, recB.manifest.rf.handedness)
        msgs{end+1} = 'RF handedness convention differs'; ok = false; %#ok<AGROW>
    end
    if recA.manifest.signal.M0 ~= recB.manifest.signal.M0
        msgs{end+1} = 'M0 normalisation differs'; ok = false; %#ok<AGROW>
    end
    if ~isequal(recA.manifest.units, recB.manifest.units)
        msgs{end+1} = 'unit declarations differ'; ok = false; %#ok<AGROW>
    end

    % ----- 3) align atoms by (T1,T2) key -----
    T1A = recA.inputs.P.T1(:); T2A = recA.inputs.P.T2(:);
    T1B = recB.inputs.P.T1(:); T2B = recB.inputs.P.T2(:);
    keyA = round([T1A T2A]*1e6);          % microsecond-resolution match key
    keyB = round([T1B T2B]*1e6);
    [tf, loc] = ismember(keyA, keyB, 'rows');

    rep.n_A       = numel(T1A);
    rep.n_B       = numel(T1B);
    rep.n_matched = nnz(tf);
    if ~all(tf)
        msgs{end+1} = sprintf('%d of %d A-atoms have no (T1,T2) match in B', ...
            nnz(~tf), numel(tf)); %#ok<AGROW>
    end

    % ----- 4) per-atom signal deviation over matched atoms -----
    rep.rmse = NaN; rep.max_abs = NaN; rep.min_correlation = NaN; rep.within_tol = false;
    rep.per_atom = struct('rmse',[],'max_abs',[],'correlation',[]);
    if rep.n_matched > 0 && isequal(shA(1), shB(1))
        MA = recA.outputs.Mxy_real + 1i*recA.outputs.Mxy_imag;
        MB = recB.outputs.Mxy_real + 1i*recB.outputs.Mxy_imag;
        iA = find(tf); iB = loc(tf);
        D  = MA(:, iA) - MB(:, iB);

        rmse_a   = sqrt(mean(abs(D).^2, 1));
        maxabs_a = max(abs(D), [], 1);
        corr_a   = zeros(1, numel(iA));
        for k = 1:numel(iA)
            a = MA(:, iA(k)); b = MB(:, iB(k));
            corr_a(k) = abs(a' * b) / (norm(a)*norm(b) + eps);
        end

        rep.per_atom.rmse        = rmse_a;
        rep.per_atom.max_abs     = maxabs_a;
        rep.per_atom.correlation = corr_a;
        rep.rmse            = sqrt(mean(abs(D(:)).^2));
        rep.max_abs         = max(abs(D(:)));
        rep.min_correlation = min(corr_a);
        rep.within_tol      = rep.max_abs <= tol;
    elseif rep.n_matched > 0
        msgs{end+1} = 'readout counts differ -> per-atom signal comparison skipped'; %#ok<AGROW>
    end

    rep.tol        = tol;
    rep.compatible = ok;
    rep.messages   = msgs;

    print_report(rep);
end

% ------------------------------------------------------------- helpers ---
function v = get_field(s, f, dflt)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = dflt; end
end

function print_report(rep)
    fprintf('\n==================  silver-record commensurability  ==================\n');
    if rep.compatible
        fprintf('structure:   COMPATIBLE\n');
    else
        fprintf('structure:   INCOMPATIBLE\n');
    end
    for k = 1:numel(rep.messages)
        fprintf('   - %s\n', rep.messages{k});
    end
    fprintf('atoms:       A=%d  B=%d  matched=%d\n', rep.n_A, rep.n_B, rep.n_matched);
    if ~isnan(rep.max_abs)
        verdict = 'FAIL'; if rep.within_tol, verdict = 'PASS'; end
        fprintf('signal:      rmse=%.3e  max|dM|=%.3e  min corr=%.6f\n', ...
                rep.rmse, rep.max_abs, rep.min_correlation);
        fprintf('tolerance:   max|dM| <= %.1e  ->  %s\n', rep.tol, verdict);
    else
        fprintf('signal:      not compared (no matched atoms or mismatched readouts)\n');
    end
    fprintf('=====================================================================\n\n');
end
