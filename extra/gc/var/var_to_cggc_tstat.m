function stat = var_to_cggc_tstat(X,V,x,p,regmode)

% Conditional GGC F-statistic
%
% The GGC test statistic is calculated for the multivariable specified by
% the vector x, conditioning on all other variables in the system specified
% by the time series X.
%
% NOTE 1: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!
%
% NOTE 2: Only the F-test works here - we cannot use the LR test because
% the least-squares model parameters in this case are not ML parameters,
% and we don't know how to estimate those!
%
% See extra/ggc_* for statistical inference. Parameters should be
% nx = size of x, nz = number of variables - size of x.

n = size(X,1);
if isempty(V)
	[~,V]  = tsdata_to_var(X,p,regmode); % full regression
else
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match time series');
end

x = x(:)'; % vectorise multivariable indices
assert(all(x >=1 & x <= n),'Some x indices out of range');

DV = diag(V);
z = 1:n; z(x) = [];
nx = length(x);
VRx  = zeros(nx,1);
for i = 1:nx
	r = [x(i) z];
	[~,VR] = tsdata_to_var(X(r,:,:),p,regmode);
	VRx(i) = VR(1,1);

end
stat = sum(VRx)/sum(DV(x)) - 1; % F-test statistic
