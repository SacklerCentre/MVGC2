function pval = mvgc_pval(stat,tstat,nx,ny,nz,p,m,N)

% Return p-values for var GC test statistics (F or likelihood-ratio chi^2)

cdf = mvgc_cdf(tstat,nx,ny,nz,p,m,N); % function handle to CDF

pval = 1-cdf(stat);
