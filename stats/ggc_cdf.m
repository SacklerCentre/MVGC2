function cdf = ggc_cdf(tstat,nx,nz,p,m,N)

% Return handle to CDF function for GGC F-statistic
%
% NOTE: Only the F-test works here - we cannot use the LR test because
% the least-squares model parameters in this case are not ML parameters,
% and we don't know how to estimate those!

assert(strcmpi(tstat,'F'),'Only F-test available for GGC');

assert(m > p,'Insufficient observations for statistical test');
M = N*(m-p); % effective number of observations
n = nx+nz;
pn = p*n;
assert(M > pn,'Insufficient observations for F-test');
d = p*nx*(nx-1);
d2 = nx*(M-pn);
sf = d2/d;  % F scaling factor
cdf = @(stat) fcdf(sf*stat,d,d2); % anonymous function handle to CDF
