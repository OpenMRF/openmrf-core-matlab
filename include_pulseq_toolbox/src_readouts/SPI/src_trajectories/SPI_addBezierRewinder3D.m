function k = SPI_addBezierRewinder3D(k, alpha, beta, N)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% Appends a cubic Bezier rewinder to a 3D k-space curve.
%
% Inputs (optional):
%   k     : Nsamples x 3, starts at 0, ends at k_end
%   alpha : Bezier control length along end tangent (default: 0.2*|k_end|)
%   beta  : Bezier control length along radial direction (default: 0.2*|k_end|)
%   N     : number of rewinder samples (default: 250)
%
% Output:
%   k     : complex vector with appended rewinder, ending at 0

k  = reshape(k, [], 3);
ke = k(end,:);
kabs = norm(ke);

% Defaults
if nargin < 4 || isempty(N)
    N = 250;
end
if nargin < 3 || isempty(beta)
    beta = 0.0;
end
if nargin < 2 || isempty(alpha)
    alpha = 0.5;
end
alpha = alpha * kabs;
beta  = beta  * kabs;

% End tangent (last non-zero segment)
tol = 1e-8 * max(vecnorm(k,2,2));
j   = size(k,1) - 1;
tk  = ke - k(j,:);
while norm(tk) <= tol && j > 1
    j  = j - 1;
    tk = ke - k(j,:);
end
if norm(tk) <= tol
    error('Cannot estimate tangent: C is (near-)constant at the end.');
end
te = tk / norm(tk);              % unit tangent at end
re = -ke / max(kabs, eps);       % unit radial direction toward center

% Bezier control points (3D)
P0 = ke;
P3 = [0 0 0];
P1 = P0 + alpha * te;            % along end tangent
P2 = P3 - beta  * re;            % shapes approach to center (radial)

% Sample Bezier curve
t = linspace(0, 1, N+1).';
u = 1 - t;

rew = (u.^3)*P0 ...
    + 3*(u.^2).*t*P1 ...
    + 3*u.*(t.^2)*P2 ...
    + (t.^3)*P3;

% Append, avoid duplicating endpoint
k = [k; rew(2:end,:)];

end
