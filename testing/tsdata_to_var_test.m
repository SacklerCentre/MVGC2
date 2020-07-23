function [A,V,E] = tsdata_to_var_test(X,p)

[n,m,N] = size(X);
assert(p < m,'too many lags');

p1 = p+1;
pn = p*n;
p1n = p1*n;

A = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')
V = NaN;
E = NaN;

X = demean(X,true); % no constant term, normalize

M = N*(m-p);

% stack lags

X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
XL = zeros(n,p,M);
for k = 1:p
	XL(:,k,:) = reshape(X(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
end
XL = reshape(XL,pn,M);         % stack lags

A = X0/XL;                     % OLS (via QR decomposition)
if isbad(A); return; end       % something went badly wrong

if nargout > 1
	E = X0-A*XL;               % residuals
	V = (E*E')/(M-1);          % residuals covariance matrix (unbiased estimator)
	if nargout > 2             % align residuals per-trial with data (lose p lags)
		E = cat(2,nan(n,p,N),reshape(E,n,m-p,N));
	end
end

A = reshape(A,n,n,p);          % so A(:,:,k) is the k-lag coefficients matrix
