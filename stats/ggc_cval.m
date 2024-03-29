function cval = ggc_cval(tstat,pcrit,nx,nz,p,m,N)

% Return critical values for GGC F-statistic
%
% NOTE: Only the F-test works here - we cannot use the LR test because
% the least-squares model parameters in this case are not ML parameters,
% and we don't know how to estimate those!

assert(strcmpi(tstat,'F'),'Only F-test available for GGC');

icdf = ggc_icdf(tstat,nx,nz,p,m,N); % function handle to inverse CDF

cval = icdf(1-pcrit);
