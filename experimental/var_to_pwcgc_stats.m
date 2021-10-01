function [stat,nullcdf,nullicdf] = var_to_pwcgc_stats(X,V,p,regmode,tstat)

% NOTE: If V is supplied, it must have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!

[n,m,N] = size(X);
if isempty(V)
	[~,V]  = tsdata_to_var(X,p,regmode);      % full regression
else
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match VAR coefficients matrix');
end

if strcmpi(tstat,'F')
	ftest = true;
elseif strcmpi(tstat,'LR')
	ftest = false;
else
	error('Unknown test statistic');
end

DV = diag(V);
if ~ftest
	LDV = log(DV);
end

stat = nan(n);
for y = 1:n
	r = [1:y-1 y+1:n]; % omit y
	[~,VR] = tsdata_to_var(X(r,:,:),p,regmode); % reduced regression
	DVR = diag(VR);
	if ftest
		stat(r,y) = DVR./DV(r) - 1;    % F-test statistic
	else
		stat(r,y) = log(DVR) - LDV(r); % likelihood-ratio test statistic
	end
end

d = p;       % Degrees of freedom
M = N*(m-p); % effective number of observations
if ftest
	d2 = M-p*n-1; % F df2
	sf = d2/d;    % F scaling factor
	nullcdf  = @(x) fcdf(sf*x,d,d2);
	nullicdf = @(x) finv(x,d,d2)/sf;
else
	sf = M;       % chi^2 scaling factor
	nullcdf  = @(x) chi2cdf(sf*x,d);
	nullicdf = @(x) chi2inv(x,d)/sf;
end

% pval = 1-nullcdf(stat);
% pcrit = mhtcorrect(pval,alpha,mhtc);
% cval = nullicdf(1-pcrit);
% sig = pval < pcrit
