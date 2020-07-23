function [F,pval] = var_to_pwcgc(A,V,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

calc_stats = nargout > 1;

if calc_stats
	assert(nargin > 3, 'must supply regression mode (same as for parameter estimates) for dual-regression stats');
	assert(~isempty(X),'must supply time series data for dual-regression stats');
	[~,m,N] = size(X);
    M  = N*(m-p);  % chi^2 scaling factor = effective number of observations
    d  = p;        % chi^2 df and F df1
    d2 = M-p*n-1;  % F df2
    K  = d2/d;     % F scaling factor
	pval.FT = nan(n);
	pval.LR = nan(n);
end

DV = diag(V);
LDV = log(DV);

F = nan(n);

FTstat = nan(n);
LRstat = nan(n);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

	[~,VR,rep] = var2riss(A,V,y,r);
    if sserror(rep,y), continue; end % check DARE report, bail out on error

    DVR = diag(VR);
    F(r,y) = log(DVR)-LDV(r);

    if calc_stats
		[~,VR] = tsdata_to_var(X(r,:,:),p,regmode); % reduced regression
		DVR = diag(VR);
        FTstat(r,y) = DVR./DV(r) - 1;       % F-test statistic
        LRstat(r,y) = log(DVR) - LDV(r);    % likelihood-ratio test statistic
    end
end

if calc_stats
    pval.FT = 1-fcdf(K*LRstat,d,d2);
    pval.LR = 1-chi2cdf(M*FTstat,d);
end
