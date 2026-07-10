function hex = silver_sha256(data)
% silver_sha256
% SHA-256 hex digest of arbitrary data, computed via the JVM (no toolboxes).
%
% ---------------------------------------------------------------- inputs ---
% data : char / string  -> hashed as its UTF-8-ish byte sequence
%        uint8 vector    -> hashed as-is
%        numeric (other) -> hashed as its raw little-endian bytes (typecast)
%
% --------------------------------------------------------------- outputs ---
% hex  : (1 x 64) char, lowercase hexadecimal SHA-256 digest
%
% Used for the record payload checksum (silver_checksum) and the pseudonymous
% machine id in the provenance block.
%
% Author: Jannik Stebani, Experimental Physics V, Wuerzburg, Germany; V0.1, 07.07.2026

    if ischar(data) || isstring(data)
        bytes = uint8(char(data));
    elseif isa(data, 'uint8')
        bytes = data(:).';
    else
        bytes = typecast(double(data(:).'), 'uint8');
    end

    md  = java.security.MessageDigest.getInstance('SHA-256');
    dig = typecast(int8(md.digest(bytes)), 'uint8');     % 32 signed->unsigned bytes
    hex = lower(reshape(dec2hex(dig, 2).', 1, []));
end
