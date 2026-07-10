function hex = silver_checksum(rec)
% silver_checksum
% Deterministic content hash of a silver record's *numeric payload*.
%
% Hashes the raw double bytes of the commensurability-relevant arrays in a
% fixed order, independent of JSON text formatting and independent of the
% descriptive metadata (manifest / provenance). Two records produced by
% different tools on different machines but carrying identical numbers
% therefore share the same payload hash -> a one-line equality check.
%
% Column-major flatten (:) is used everywhere so the hash is invariant to the
% row/column orientation that jsondecode may hand back for 1-D arrays.
%
% ---------------------------------------------------------------- inputs ---
% rec : silver record struct (see make_silver_record). Requires the numeric
%       payload fields under rec.inputs and rec.outputs.
%
% --------------------------------------------------------------- outputs ---
% hex : (1 x 64) char lowercase SHA-256 of the concatenated payload doubles.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    in  = rec.inputs;
    out = rec.outputs;
    parts = { in.FAs_real(:).', in.FAs_imag(:).', ...
              in.TRs(:).',      in.TE(:).', ...
              in.P.T1(:).',     in.P.T2(:).', ...
              out.Mxy_real(:).', out.Mxy_imag(:).' };
    v   = double([parts{:}]);
    hex = silver_sha256(typecast(v, 'uint8'));
end
