function stat = var_to_mvgc_tstat(X,V,x,y,p,regmode,tstat)

% Return multivariate Granger causality test statistic (F or likelihood ratio)
%
% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!

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

x = x(:)'; % vectorise target variable indices
y = y(:)'; % vectorise source variable indices

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'Some x indices out of range');
assert(all(y >=1 & y <= n),'Some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of conditioning variables (i.e., all other variablesin model)
r = [x z];               % indices of variables in reduced model (omit source variables)

[~,VR]  = tsdata_to_var(X(r,:,:),p,regmode);  % reduced regression

xr = 1:nx; % indices of target in reduced model

if ftest
	stat = trace(VR(xr,xr))/trace(V(x,x)) - 1;  % F-test statistic
else
	stat = logdet(VR(xr,xr)) - logdet(V(x,x));  % likelihood-ratio test statistic
end
