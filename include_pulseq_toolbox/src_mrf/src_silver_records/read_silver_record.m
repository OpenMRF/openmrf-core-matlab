function rec = read_silver_record(path, verify)
% read_silver_record
% Load a silver record from JSON, reconstruct its complex arrays and (by
% default) verify the payload checksum.
%
% ---------------------------------------------------------------- inputs ---
% path   : path to a .json silver record written by write_silver_record.
% verify : (optional) logical, re-check integrity.payload_sha256.  default true
%
% --------------------------------------------------------------- outputs ---
% rec : silver record struct. In addition to the on-disk fields it carries two
%       convenience reconstructions:
%         rec.Mxy : (NR x N_dict) complex dictionary
%         rec.FAs : (NR x 1)      complex flip-angle pattern
%       Mxy is reshaped to the manifest-declared [NR N_dict] shape, so the
%       result is robust to jsondecode's 1-D orientation quirks.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    if nargin < 2 || isempty(verify), verify = true; end

    txt = fileread(path);
    rec = jsondecode(txt);

    % ----- reshape / orient the numeric payload to the declared layout -----
    shp = double(rec.manifest.dims.shape(:).');
    rec.outputs.Mxy_real = reshape(rec.outputs.Mxy_real, shp);
    rec.outputs.Mxy_imag = reshape(rec.outputs.Mxy_imag, shp);

    rec.inputs.FAs_real = rec.inputs.FAs_real(:);
    rec.inputs.FAs_imag = rec.inputs.FAs_imag(:);
    rec.inputs.TRs      = rec.inputs.TRs(:);
    rec.inputs.TE       = rec.inputs.TE(:);
    rec.inputs.P.T1     = rec.inputs.P.T1(:);
    rec.inputs.P.T2     = rec.inputs.P.T2(:);

    % ----- integrity check -----
    if verify && isfield(rec, 'integrity') && isfield(rec.integrity, 'payload_sha256')
        got = silver_checksum(rec);
        if ~strcmp(got, rec.integrity.payload_sha256)
            warning('read_silver_record:checksum', ...
                ['payload checksum mismatch (file corrupted or hand-edited).\n' ...
                 '  stored: %s\n  actual: %s'], rec.integrity.payload_sha256, got);
        end
    end

    % ----- convenience complex reconstructions -----
    rec.Mxy = rec.outputs.Mxy_real + 1i*rec.outputs.Mxy_imag;
    rec.FAs = rec.inputs.FAs_real  + 1i*rec.inputs.FAs_imag;
end
