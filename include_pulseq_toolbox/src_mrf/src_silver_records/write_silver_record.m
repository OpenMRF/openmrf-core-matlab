function outpath = write_silver_record(rec, outpath)
% write_silver_record
% Serialize a silver record to a single self-contained JSON file.
%
% Complex arrays are stored split into *_real / *_imag (see the manifest); the
% file is fully self-describing and language-agnostic. Reload with
% read_silver_record.
%
% ---------------------------------------------------------------- inputs ---
% rec     : silver record struct from make_silver_record.
% outpath : target path. A missing extension gets '.json'. Parent
%           directories are created if needed.
%
% --------------------------------------------------------------- outputs ---
% outpath : the (possibly extension-completed) path actually written.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    if nargin < 2 || isempty(outpath)
        error('write_silver_record:path', 'an output path is required.');
    end

    [d, ~, ext] = fileparts(outpath);
    if isempty(ext), outpath = [outpath '.json']; end
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

    try
        txt = jsonencode(rec, 'PrettyPrint', true);   % R2021a+
    catch
        txt = jsonencode(rec);                         % older MATLAB
    end

    fid = fopen(outpath, 'w');
    if fid < 0, error('write_silver_record:open', 'cannot open %s for writing.', outpath); end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid, txt, 'char');

    fprintf('silver record written: %s  (%d readouts x %d tissues)\n', ...
            outpath, rec.manifest.dims.shape(1), rec.manifest.dims.shape(2));
end
