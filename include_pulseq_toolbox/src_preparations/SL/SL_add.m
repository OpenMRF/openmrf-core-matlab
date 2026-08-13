% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% temporary SL objects
exc1     = SL.SL_objs(loop_SL).EXC1;
exc2     = SL.SL_objs(loop_SL).EXC2;
sl1      = SL.SL_objs(loop_SL).SL1;
sl2      = SL.SL_objs(loop_SL).SL2;
sl3      = SL.SL_objs(loop_SL).SL3;
rfc1     = SL.SL_objs(loop_SL).RFC1;
rfc2     = SL.SL_objs(loop_SL).RFC2;
d1       = SL.SL_objs(loop_SL).d1;
d2       = SL.SL_objs(loop_SL).d2;
gx_crush = SL.SL_objs(loop_SL).gx_crush;
gy_crush = SL.SL_objs(loop_SL).gy_crush;
gz_crush = SL.SL_objs(loop_SL).gz_crush;

%% add excitation pulse: tip down
seq.addBlock(exc1);

%% add spin-lock cluster
seq.addBlock(d1);

if strcmp(SL.seq_type(loop_SL), 'SSL')
    seq.addBlock(sl1);    seq.addBlock(d1);

elseif strcmp(SL.seq_type(loop_SL), 'hSSL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);

elseif strcmp(SL.seq_type(loop_SL), 'qSSL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);
    
elseif strcmp(SL.seq_type(loop_SL), 'RESL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);

elseif strcmp(SL.seq_type(loop_SL), 'CSL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(rfc1);   seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);

elseif strcmp(SL.seq_type(loop_SL), 'PSCSL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);   
    seq.addBlock(rfc1);   seq.addBlock(d1);
    seq.addBlock(sl1);    seq.addBlock(d1);   
    seq.addBlock(sl2);    seq.addBlock(d1);  

elseif strcmp(SL.seq_type(loop_SL), 'TBSL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(rfc1);   seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);
    seq.addBlock(rfc2);   seq.addBlock(d1);
    seq.addBlock(sl1);    seq.addBlock(d1);

elseif strcmp(SL.seq_type(loop_SL), 'BSL')
    seq.addBlock(sl1);    seq.addBlock(d1);
    seq.addBlock(rfc1);   seq.addBlock(d1);
    seq.addBlock(sl2);    seq.addBlock(d1);
    seq.addBlock(rfc2);   seq.addBlock(d1);
    seq.addBlock(sl3);    seq.addBlock(d1);

end

%% add excitation pulse: tip up
if ~isempty(exc2)
    seq.addBlock(exc2);
end

%% add crusher gradients
seq.addBlock(gx_crush, gy_crush, gz_crush);
seq.addBlock(d2);                

%%
clear exc1 exc2 sl1 sl2 sl3 rfc1 rfc2 d1 d2 gx_crush gy_crush gz_crush;