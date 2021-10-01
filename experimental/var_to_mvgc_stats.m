function [stat,nullcdf,nullicdf] = var_to_mvgc_stats(X,V,p,x,y,regmode,tstat)

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

x = x(:)'; % vectorise target variable indices
y = y(:)'; % vectorise source variable indices

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'Some x indices out of range');
assert(all(y >=1 & y <= n),'Some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of conditioning variables (i.e., all other variablesin model)
r = [x z];               % indices of variables in reduced model (omit source variables)

[~,VR]  = tsdata_to_var(X(r,:,:),p,regmode);  % reduced regression

nx = length(x);
ny = length(y);
xr = 1:nx; % indices of target in reduced model

d = p*nx*ny; % Degrees of freedom
M = N*(m-p); % effective number of observations
if ftest
	d2 = nx*(M-p*n)-1; % F df2
	sf = d2/d;         % F scaling factor
	stat = trace(VR(xr,xr))/trace(V(x,x)) - 1;  % F-test statistic
	nullcdf  = @(x) fcdf(sf*x,d,d2);
	nullicdf = @(x) finv(x,d,d2)/sf;
else
	sf = M;            % chi^2 scaling factor = effective number of observations
	stat = logdet(VR(xr,xr)) - logdet(V(x,x));  % likelihood-ratio test statistic
	nullcdf  = @(x) chi2cdf(sf*x,d);
	nullicdf = @(x) chi2inv(x,d)/sf;
end

% pval = 1-nullcdf(stat);
% pcrit = mhtcorrect(pval,alpha,mhtc);
% cval = nullicdf(1-pcrit);
% sig = pval < pcrit
