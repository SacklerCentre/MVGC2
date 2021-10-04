function [pval,stats] = pwcgc_var_stats(X,V,p,regmode,tstat)

% Return pairwise-continuous Granger causality test statistics
% (F or likelihood ratio) and p-values.

% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!

[n,m,N] = size(X);
if isempty(V)
	[~,V]  = tsdata_to_var(X,p,regmode); % full regression
else
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match time series');
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

stats = nan(n);
for y = 1:n
	r = [1:y-1 y+1:n]; % omit y
	[~,VR] = tsdata_to_var(X(r,:,:),p,regmode); % reduced regression
	DVR = diag(VR);
	if ftest
		stats(r,y) = DVR./DV(r) - 1;    % F-test statistic
	else
		stats(r,y) = log(DVR) - LDV(r); % likelihood-ratio test statistic
	end
end

d = p;       % Degrees of freedom
M = N*(m-p); % effective number of observations
if ftest
	d2 = M-p*n-1; % F df2
	sf = d2/d;    % F scaling factor
	pval = 1-fcdf(sf*stats,d,d2);
else
	sf = M;       % chi^2 scaling factor
	pval = 1-chi2cdf(sf*stats,d);
end
