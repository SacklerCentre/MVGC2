function [F,pval] = var_to_mvgc(A,V,x,y,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

calc_stats = nargout > 1;

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];               % indices of reduced variables

nx = length(x);
ny = length(y);
xr = 1:nx;               % index of x in reduced variables

if calc_stats
	assert(nargin > 5, 'must supply regression mode (same as for parameter estimates) for dual-regression stats');
	assert(~isempty(X),'must supply time series data for dual-regression stats');
	[~,m,N] = size(X);
	M  = N*(m-p);      % chi^2 scaling factor = effective number of observations
	d  = p*nx*ny;      % chi^2 df and F df1
	d2 = nx*(M-p*n)-1; % F df2
	K  = d2/d;         % F scaling factor
end

F = NaN;

[~,VR,rep] = var2riss(A,V,y,r);
if sserror(rep), return; end % check DARE report, bail out on error

F = logdet(VR(xr,xr)) - logdet(V(x,x));

if calc_stats
	[~,VR]  = tsdata_to_var(X(r,:,:),p,regmode);  % reduced regression
	FTstat  = trace(VR(xr,xr))/trace(V(x,x)) - 1; % F-test statistic
	pval.FT = 1-fcdf(K*FTstat,d,d2);
	LRstat  = logdet(VR(xr,xr)) - logdet(V(x,x)); % likelihood-ratio test statistic
	pval.LR = 1-chi2cdf(M*LRstat,d);
end
