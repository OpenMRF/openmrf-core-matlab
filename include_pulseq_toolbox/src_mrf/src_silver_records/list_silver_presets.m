function T = list_silver_presets()
% list_silver_presets
% Enumerate the available record-flavour presets by reading their JSON
% definitions in ./presets/. Reads metadata only -- runs no simulation.
%
% --------------------------------------------------------------- outputs ---
% T : table with columns name, alias, mode, title. A human-readable listing
%     is also printed to the console.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    d     = fullfile(fileparts(mfilename('fullpath')), 'presets');
    files = dir(fullfile(d, '*.json'));

    rows = struct('name', {}, 'alias', {}, 'mode', {}, 'title', {});
    for k = 1:numel(files)
        c = jsondecode(fileread(fullfile(d, files(k).name)));
        rows(end+1) = struct( ...
            'name',  getf(c, 'name',  ''), ...
            'alias', getf(c, 'alias', ''), ...
            'mode',  getf(getf(c, 'tissues', struct()), 'mode', ''), ...
            'title', getf(c, 'title', '')); %#ok<AGROW>
    end

    fprintf('\n available silver-record presets (%d):\n', numel(rows));
    for k = 1:numel(rows)
        a = rows(k).alias;
        if ~isempty(a), a = sprintf('  (alias: %s)', a); end
        fprintf('   %-22s%s\n      %s\n', rows(k).name, a, rows(k).title);
    end
    fprintf('\n');

    if isempty(rows)
        T = table();
    else
        T = struct2table(rows);
    end
end

function v = getf(s, f, dflt)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
