function [F,pval] = var_to_pwcgc(A,V,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

DV = diag(V);
LDV = log(DV);

F = nan(n);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
	[~,VR,rep] = var2riss(A,V,y,r);
    if sserror(rep,y), continue; end % check DARE report, bail out on error
    F(r,y) = log(diag(VR))-LDV(r);
end

if nargout > 1 % calculate stats
	assert(nargin > 3, 'Must supply regression mode for stats (same mode as used for VAR model estimate)');
	assert(~isempty(X),'Must supply time-series data for stats');
	[n1,m,N] = size(X);
	assert(n1 == n,    'Time series does not match VAR coefficients matrix');
    M  = N*(m-p);  % chi^2 scaling factor = effective number of observations
    d  = p;        % chi^2 df and F df1
    d2 = M-p*n-1;  % F df2
    K  = d2/d;     % F scaling factor
	FTstat = nan(n);
	LRstat = nan(n);
	for y = 1:n
		r = [1:y-1 y+1:n]; % omit y
		[~,VR] = tsdata_to_var(X(r,:,:),p,regmode); % reduced regression
		DVR = diag(VR);
        FTstat(r,y) = DVR./DV(r) - 1;    % F-test statistic
        LRstat(r,y) = log(DVR) - LDV(r); % likelihood-ratio test statistic
    end
	pval.FT = nan(n);
	pval.LR = nan(n);
    pval.FT = 1-fcdf(K*LRstat,d,d2);
    pval.LR = 1-chi2cdf(M*FTstat,d);
end
