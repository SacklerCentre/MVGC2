function cval = mvgc_cval(pcrit,tstat,nx,ny,nz,p,m,N)

% Return critical values for var GC test statistics (F or likelihood-ratio chi^2)

icdf = mvgc_icdf(tstat,nx,ny,nz,p,m,N); % function handle to inverse CDF

cval = icdf(1-pcrit);
