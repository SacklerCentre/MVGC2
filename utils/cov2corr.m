% Convert covariance matrix to correlation matrix
%
% Result guaranteed positive-definite if covariance matrix is
% positive-definite (it should be! - check p == 0), otherwise
% a result is still returned, but should be treated as suspect.

function [R,p,pval] = cov2corr(V,N)

[n,n1] = size(V);
assert(ismatrix(V) && n1 == n,'Covariance matrix must be square');

R  = NaN(n);
pval = NaN(n);

[R,p] = chol(V); % right Cholesky factor
if p > 0 % fall back on "safe" method (doesn't ensure positive-definite result)
	d = 1./sqrt(diag(V));
	R = bsxfun(@times,bsxfun(@times,d,V),d');
else      % ensures positive-definite result
	R = bsxfun(@rdivide,R,sqrt(sum(R.*R)));
	R = R'*R;
end

if nargout > 2 % calculate (approximate) p-values using Fisher z-transformation
	assert(nargin > 1 && ~isempty(N),'Need sample size for p-values');
	z = sqrt(N-3)*atanh(R);           % z-score for Fisher transformation
	pval = 2*(1-normcdf(abs(z)));     % 2-tailed test
end
