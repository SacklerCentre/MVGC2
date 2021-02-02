function [F,pval] = var_to_mvgc(A,V,x,y,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match VAR coefficients matrix');

x = x(:)'; % vectorise target variable indices
y = y(:)'; % vectorise source variable indices

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'Some x indices out of range');
assert(all(y >=1 & y <= n),'Some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of conditioning variables (i.e., all other variablesin model)
r = [x z];               % indices of variables in reduced model (omit source variables)

nx = length(x);
ny = length(y);
xr = 1:nx;               % indices of target in reduced model

F = NaN;
[~,VR,rep] = var2riss(A,V,y,r);
if sserror(rep), return; end % check DARE report, bail out on error
F = logdet(VR(xr,xr)) - logdet(V(x,x));

if nargout > 1 % calculate stats
	assert(nargin > 5, 'Must supply regression mode for stats (same mode as used for VAR model estimate)');
	assert(~isempty(X),'Must supply time-series data for stats');
	[n1,m,N] = size(X);
	assert(n1 == n,    'Time series does not match VAR coefficients matrix');
	M  = N*(m-p);      % chi^2 scaling factor = effective number of observations
	d  = p*nx*ny;      % chi^2 df and F df1
	d2 = nx*(M-p*n)-1; % F df2
	K  = d2/d;         % F scaling factor
	[~,VR]  = tsdata_to_var(X(r,:,:),p,regmode);  % reduced regression
	FTstat  = trace(VR(xr,xr))/trace(V(x,x)) - 1; % F-test statistic
	LRstat  = logdet(VR(xr,xr)) - logdet(V(x,x)); % likelihood-ratio test statistic
	pval.FT = 1-fcdf(K*FTstat,d,d2);
	pval.LR = 1-chi2cdf(M*LRstat,d);
end
