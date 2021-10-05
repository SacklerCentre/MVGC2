function pval = mvgc_pval(stat,tstat,nx,ny,nz,p,m,N)

% Return p-values for var GC test statistics (F or likelihood-ratio chi^2)

if strcmpi(tstat,'F')
	ftest = true;
elseif strcmpi(tstat,'LR')
	ftest = false;
else
	error('Unknown test statistic');
end

d = p*nx*ny; % Degrees of freedom
M = N*(m-p); % effective number of observations
if ftest
	n = nx+ny+nz;
	d2 = nx*(M-p*n)-1; % F df2
	sf = d2/d;         % F scaling factor
	pval = 1-fcdf(sf*stat,d,d2);
else
	sf = M;            % chi^2 scaling factor
	pval = 1-chi2cdf(sf*stat,d);
end
