function P = mvgc_F_cdf(x,p,m,N,nx,ny,nz)

% For use with Granger form (F distribution), although strictly speaking only
% really works with dual regression estimation, hence not that useful. Also,
% only works for the null distribution.

xx = x; xx(isnan(x)) = [];
assert(all(xx >= 0), 'MVGC sample values must be non-negative');

if isempty(N), N = 1; end % single trial

if nargin < 7 || isempty(nz), nz = 0; end % unconditional

m  = N*(m-p);             % effective number of observations
d1 = p*nx*ny;             % F df1
d2 = nx*(m-p*(nx+ny+nz)); % F df2
k  = d2/d1;               % F scaling factor
P  = fcdf(k*x,d1,d2);     % theoretical F cdf
