function pval = ggc_pval(stat,tstat,nx,nz,p,m,N)

% Return p-values for GGC F-statistic
%
% NOTE: Only the F-test works here - we cannot use the LR test because
% the least-squares model parameters in this case are not ML parameters,
% and we don't know how to estimate those!

assert(strcmpi(tstat,'F'),'Only F-test available for GGC');

cdf = ggc_cdf(tstat,nx,nz,p,m,N); % function handle to CDF

pval = 1-cdf(stat);
