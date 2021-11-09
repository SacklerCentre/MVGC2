function cdf = mvgc_cdf(tstat,nx,ny,nz,p,m,N)

% Return handle to CDF function for GC test statistics (F or likelihood-ratio chi^2)

if strcmpi(tstat,'F')
	ftest = true;
elseif strcmpi(tstat,'LR')
	ftest = false;
else
	error('Unknown test statistic');
end

assert(m > p,'Insufficient observations for statistical tests');
d = p*nx*ny; % degrees of freedom
M = N*(m-p); % effective number of observations
if ftest
	n = nx+ny+nz;
	pn = p*n;
	assert(M > pn,'Insufficient observations for F-test');
	d2 = nx*(M-pn); % F df2
	sf = d2/d;      % F scale factor
	cdf = @(stat) fcdf(sf*stat,d,d2); % anonymous function handle to CDF
else
	sf = M;         % chi^2 scale factor
	cdf = @(stat) chi2cdf(sf*stat,d); % anonymous function handle to CDF
end
