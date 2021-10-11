function bias = mvgc_bias(tstat,nx,ny,nz,p,m,N)

% Return approximate bias (scaled null mean) of MVGC test statistics (F or likelihood-ratio chi^2)

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
	bias = d/(d2-2);   % scaled mean of F distribution
else
	bias = d/M;        % scaled mean of chi^2 distribution
end
