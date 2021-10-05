function cval = mvgc_cval(pcrit,tstat,nx,ny,nz,p,m,N)

% Return critical values for var GC test statistics (F or likelihood-ratio chi^2)

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
	cval = finv(1-pcrit,d,d2)/sf;
else
	sf = M;            % chi^2 scaling factor
	cval = chi2inv(1-pcrit,d)/sf;
end
