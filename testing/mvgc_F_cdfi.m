function x = mvgc_F_cdfi(P,p,m,N,nx,ny,nz)

% For use with Granger form (F distribution), although strictly speaking only
% really works with dual regression estimation, hence not that useful. Also,
% only works for the null distribution.

PP = P; PP(isnan(P)) = [];
assert(all(PP >= 0) && all(PP <= 1), 'probabilities must lie between 0 and 1');

if isempty(N), N = 1; end % single trial

if nargin < 7 || isempty(nz), nz = 0; end % unconditional

m  = N*(m-p);             % effective number of observations
d1 = p*nx*ny;             % F df1
d2 = nx*(m-p*(nx+ny+nz)); % F df2
k  = d2/d1;               % F scaling factor
x  = finv(P,d1,d2)/k;     % theoretical F inverse cdf
