function stat = var_to_pwcgc_tstat(X,V,p,regmode,tstat)

% Return pairwise-continuous Granger causality test statistics (F or likelihood ratio) and p-values.

% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!
%
% See stats/mvgc_* for statistical inference. Parameters should be
% nx = 1, ny = 1, nz = total number of variables - 2.

if strcmpi(tstat,'F')
	ftest = true;
elseif strcmpi(tstat,'LR')
	ftest = false;
else
	error('Unknown test statistic');
end

n = size(X,1);
if isempty(V)
	[~,V]  = tsdata_to_var(X,p,regmode); % full regression
else
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match time series');
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
