function [R,p,pval] = cov2corr(V,N)

% Calculate pairwise (Pearson) correlation matrix from covariance
% matrix. P-values are calculated for a 2-tailed test.
%
% The algorithm attempts to ensure a positive-definite result; if
% the covariance matrix is not positive-definite, we fall back on
% a "safe" method.
%
% V - sample covariance matrix (unbiased estimator)
% N - sample size (number of data observations)

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
