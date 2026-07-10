function [P, labels] = default_silver_tissues(names)
% default_silver_tissues
% A small, curated set of nominal (T1,T2) tissue relaxation times for building
% reproducible silver records. Values are representative ~3T literature figures
% intended as stable validation anchors -- NOT clinical ground truth.
%
% ---------------------------------------------------------------- inputs ---
% names : (optional) char or cellstr of tissue names to select (in order).
%         Omit / empty -> return the full table.
%
% --------------------------------------------------------------- outputs ---
% P      : struct with P.T1, P.T2 (N x 1) [s], ready for make_silver_record.
% labels : {N x 1} cellstr of the selected tissue names.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    % name                 T1 [s]   T2 [s]   (nominal, ~3T)
    tbl = {
        'white_matter',    0.830,   0.080
        'gray_matter',     1.330,   0.110
        'csf',             4.000,   2.000
        'skeletal_muscle', 1.230,   0.037
        'fat',             0.370,   0.130
        'blood_arterial',  1.660,   0.150
        'myocardium',      1.180,   0.046
        'liver',           0.810,   0.042
        'cartilage',       1.240,   0.030
        'kidney_medulla',  1.550,   0.081
    };
    all_names = tbl(:,1);

    if nargin < 1 || isempty(names)
        sel = (1:size(tbl,1))';
    else
        if ischar(names), names = cellstr(names); end
        sel = zeros(numel(names), 1);
        for k = 1:numel(names)
            j = find(strcmpi(all_names, names{k}), 1);
            if isempty(j)
                error('default_silver_tissues:unknown', 'unknown tissue: %s', names{k});
            end
            sel(k) = j;
        end
    end

    labels = tbl(sel, 1);
    P.T1   = cell2mat(tbl(sel, 2));
    P.T2   = cell2mat(tbl(sel, 3));
    P.T1   = P.T1(:);
    P.T2   = P.T2(:);
end
