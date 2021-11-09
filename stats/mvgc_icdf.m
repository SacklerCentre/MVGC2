function icdf = mvgc_icdf(tstat,nx,ny,nz,p,m,N)

% Return handle to inverse CDF function for GC test statistics (F or likelihood-ratio chi^2)

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
	icdf = @(prob) finv(prob,d,d2)/sf; % anonymous function handle to inverse CDF
else
	sf = M;         % chi^2 scale factor
	icdf = @(prob) chi2inv(prob,d)/sf; % anonymous function handle to inverse CDF
end
