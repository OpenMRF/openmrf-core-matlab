function M_dict = MRF_sim_EPG(SIM, P, comp_energy)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% this function simulates a dictionary M_dict for different parameters P
% the simulation is based on the epg formalism
% rf pulses are simulated as instantaneous events
% relaxation losses during adiabatic pulses are corrected
% the simulation models:
%   T1/T2 relaxation
%   T2' relaxation
%   T1p relaxation (neglecting T2p)
%   B0 and B1+ deviations
%   apparent diffusion
% approximations:
%   unit gradients
%   ideal slice profiles
%   adiabatic losses approximated via look-up-tables

% ----- input -----
% SIM
%   .ID :  (N_steps x 1) ID for Bloch operator [1...10]
%   .RF :  (N_steps x 1) complex RF waveform [rad/s]
%   .GZ :  (N_steps x 1) z gradient moments [rad/m] or waveforms [rad/s/m]
%   .DB :  (N_steps x 1) b-value of the diffusion encoding tensor [s/mm^2]
%   .DT :  (N_steps x 1) dynamic time steps [s]
%   .PHI:  (N_adc x 1)   phase of receiver [rad]
% P
%   .T1:  (N_dict x 1) T1 relaxation time [s]
%   .T2:  (N_dict x 1) T2 relaxation time [s]
%   .T1p: (N_dict x 1) T1p relaxation time [s]
%   .m1p: (N_dict x 1) linear T1p dispersion [ms/kHz]
%   .T1p_adia: (N_dict x 1) adiabatic T1p relaxation time [s]
%   .adc  (N_dict x 1) apparent diffusion coefficient (mm^2/s)
%   .db1: (N_dict x 1) B1+ field deviation
%   .dw0: (N_dict x 1) B0 field deviation [rad/s]
%   .TM:  (N x 3x3)    pre-computed transition matrix for adiabatic pulses
% comp_energy: 0 -> full dict; >0 -> compressed dict (e.g. 99.99%)

% ----- output -----
% M_dict : (N_adc x N_dict) complex signal

% load pre-compressed dictionary
if nargin<3
    comp_energy = 0; % no compression
end
if comp_energy>0
    [dict_comp, dict_hash] = MRF_find_comp_dict(comp_energy, SIM, P);
    if ~isempty(dict_comp)
        M_dict = dict_comp;
        fprintf(' \n');
        fprintf(['loading compressed dictionary from: ' dict_hash ' \n']);
        fprintf(' \n');
        return;
    end
end

% default parameters
N_dict = numel(P.T1);
if isfield(P, 'T1')
    T1 = P.T1;
else
    error('define T1!');
end
if isfield(P, 'T2')
    T2 = P.T2;
else
    error('define T2!');
end
if isfield(P, 'db1')
    db1 = P.db1;
else
    db1 = ones(N_dict, 1);
end
if isfield(P, 'dw0')
    dw0 = P.dw0;
else
    dw0 = zeros(N_dict, 1);
end
if isfield(P, 'T1p')
    T1p = P.T1p;
else
    T1p = T2;
end
if isfield(P, 'm1p')
    m1p = P.m1p;
else
    m1p = zeros(N_dict, 1);
end
if isfield(P, 'T1p_adia')
    T1p_adia = P.T1p_adia;
else
    T1p_adia = T1p;
end
if isfield(P, 'ADC')
    ADC = P.ADC;
else
    ADC = zeros(N_dict, 1);
end
clear P;

% operator lists
ID  = SIM.ID;
RF  = SIM.RF;
GZ  = SIM.GZ;
DB  = SIM.DB;
DT  = SIM.DT;
PHI = SIM.PHI;
TM  = SIM.TM;
clear SIM;

% init
N_steps = numel(ID);
N_adc   = sum(ID==0);
N_q     = ceil(sum(abs(GZ))/2);
M_dict  = zeros(N_adc, N_dict);

tic;
fprintf('\n');
fprintf(['calculating dictionary: ' num2str(N_dict) ' atoms via EPG \n']);

% calculate dictionary
parfor j = 1:N_dict

    % init magnetization, Mz = 1
    Q      = zeros(3, N_q);
    Q(3,1) = 1;
    M_sim  = zeros(N_adc, 1);

    % B1+ deviation
    RF_ = RF .* (db1(j)*(ID~=3) + (ID==3)); % apply for all RF events, except adiabatic pulses

    % counter for ADCs
    n_adc = 1;
    
    % forward simulation using different Bloch operators
    for k = 1:N_steps
    switch(ID(k))
    
        case 0 % adc
            M_sim(n_adc) = Q(1,1);
            M_sim(n_adc) = M_sim(n_adc) * exp(-1i * PHI(n_adc));
            n_adc = n_adc + 1;               

        case 1 % free relaxation & dephasing
            Q = sim_free_relax(Q, T1(j), T2(j), DT(k));
            Q = sim_dephasing(Q, DT(k)*dw0(j));

        case {2, 5, 6} % rf excitation
            Q = sim_rf(Q, abs(RF_(k)), angle(RF_(k)));

        case 3 % adiabatic rf excitation
            Q = squeeze(TM(j,abs(RF_(k)),:,:)) * Q;

        case 4 % spin-lock
            Q = sim_sl(Q, DT(k), RF_(k), T1p(j), m1p(j));

        case {8, 9} % perfect crusher
            Q = [0 0 0; 0 0 0; 0 0 1] * Q;

        case 10 % z gradient
            Q = sim_spoiler(Q, GZ(k));

        case 11 % diffusion
            Q = sim_diffusion(Q, DB(k), ADC(j));

        case 12 % adiabatic spin-lock
            Q = sim_rf(Q, pi, 0);
            Q = exp(-DT(k)/T1p_adia(j)) * Q;
    end
    end

    M_dict(:,j) = M_sim;

end

temp_t = toc;
fprintf(['   ' num2str(temp_t, '%.1f') 's ... complete! '  num2str(temp_t/N_dict*1e3, '%.1f')  ' ms/atom \n']);

% save pre-compressed dictionary
if comp_energy>0
    [~, S ,V] = svd(M_dict.', 'econ');
    svals     = diag(S)/S(1);
    svals     = cumsum(svals.^2./sum(svals.^2));
    NPCs      = find(svals > comp_energy, 1, 'first');
    dict_phi  = single(V(:,1:NPCs));
    dict_svd  = (M_dict.'*dict_phi).';
    dict_comp.phi  = dict_phi;
    dict_comp.svd  = dict_svd;
    dict_comp.NPCs = NPCs;
    save(dict_hash, 'dict_comp', '-v7.3');
    M_dict = dict_comp;
end

end

%% ------------------- EPG operator: free relaxation -------------------
function Q = sim_free_relax(Q, T1, T2, dt)
 
% Weigel M. Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple.
% Journal of Magnetic Resonance Imaging. 2014;41(2):266-295. doi:10.1002/jmri.24619
% Eq. 24

% input:
% Q matrix old, 3xN
% T1, T2 relaxation times [s]
% dt, duration of free relaxation [s]

% output
%Q matrix new, 3xN

Q = [  exp(-dt/T2)  0            0;
       0            exp(-dt/T2)  0;
       0            0            exp(-dt/T1)
    ] * Q;

Q(:,1) = Q(:,1) + [0; 0; 1-exp(-dt/T1)];

end

%% ------------------- EPG operator: free dephasing -------------------
function Q = sim_dephasing(Q, dphi)

% Weigel M. Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple.
% Journal of Magnetic Resonance Imaging. 2014;41(2):266-295. doi:10.1002/jmri.24619
% Eq. 11

% input:
% Q matrix old, 3xN
% dphi, depahsing angle [rad], e.g. dphi = dw0 * tau

% output
% Q matrix new, 3xN

Q = [  exp(-1i*dphi)   0               0;
       0               exp(1i*dphi)    0;
       0               0               1
    ] * Q;

end

%% ------------------- EPG operator: rf pulses -------------------
function Q = sim_rf(Q, alpha, phi)

% Weigel M. Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple.
% Journal of Magnetic Resonance Imaging. 2014;41(2):266-295. doi:10.1002/jmri.24619
% Eq. 18

% input:
% Q matrix old, 3xN
% alpha, flip angle [rad]
% phi, roation axis [rad]

% output
%Q matrix new, 3xN

alpha = -alpha; % M. Gram: right-hand -> left-hand rotation!

Q = [  cos(alpha/2)^2                       exp(2*1i*phi) * sin(alpha/2)^2     -1i * exp( 1i*phi) * sin(alpha);
       exp(-2*1i*phi) * sin(alpha/2)^2      cos(alpha/2)^2                      1i * exp(-1i*phi) * sin(alpha);
       -1i/2 * exp(-1i*phi) * sin(alpha)    1i/2 * exp(1i*phi) * sin(alpha)     cos(alpha)
    ] * Q;

end

%% ------------------- EPG operator: gradient spoiler -------------------
function Q = sim_spoiler(Q, nTwists)

N = size(Q, 2); % number of q-states

% positive gradient: dephasing
if nTwists>0
    for n=1:abs(nTwists)        
        for j = 1:N-1
            Q(1, N - j + 1) = Q(1, N - j);
            Q(2, j)         = Q(2, j + 1);
        end
        Q(2, N) = 0;
        Q(1, 1) = conj(Q(2, 1));
    end
end

% negative gradient: rephasing
if nTwists<0
    for n=1:abs(nTwists)
        Q(2, 1) = Q(1, 1);
        for j = 1:N-1
            Q(1, j)         = Q(1, j + 1);
            Q(2, N - j + 1) = Q(2, N - j);
        end
        Q(1, N) = 0;
        Q(2, 1) = conj(Q(1,1));
    end
end

end

%% ------------------- EPG operator: diffusion -------------------
function Q = sim_diffusion(Q, B, ADC)
    B = [ exp(-B*ADC), 0, 0,;
          0, exp(-B*ADC), 0;
          0, 0,           1 ];
    Q = B * Q;
end

%% ------------------- EPG operator: spin-lock -------------------
function Q = sim_sl(Q, tSL, wSL, T1p, m1p)

% M. Gram: epg simulation for spin-locking
% ... it's a problem
% the Q state formalism does not easily allow to simulate spin-lock conditions
% consider that the relaxation operator for Q states can not be adapted for T1p
% and T2p since the spin-lock forces new axis for relaxation processes
% ... as an approximation, I will subdivide the effect of spin-locking in
% I)   relaxation: we neglect T2p effects and attenuate all Q states with T1p decay
% II)  rotation:   the effective flip angle of the spin-lock is simulated as an rf operator

% assume linear T1p dispersion
T1p = T1p + m1p * abs(wSL) * 1e-6 / 2 / pi; % m1p -> [ms/kHz]

% simulate T1p decay
Q = Q * exp(-tSL/T1p);

% simulate spin-lock-rotation
Q = sim_rf(Q, abs(wSL)*tSL, angle(wSL));

end