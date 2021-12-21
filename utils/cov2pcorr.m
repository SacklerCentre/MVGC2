function [R,p,pval] = cov2pcorr(V,N)

% Calculate pairwise partial (Pearson) correlation coefficients from
% covariance matrix using inversion method. P-values are calculated
% for a 2-tailed test.
%
% V - sample covariance matrix (unbiased estimator)
% N - sample size (number of data observations)
%
% NOTE: The pairwise-conditional mutual information under Gaussian
% assumptions is -log(1-R.*conj(R)), with an asymptotic chi^2(1)
% distribution.

[n,n1] = size(V);
assert(ismatrix(V) && n1 == n,'Covariance matrix must be square');

R  = NaN(n);
pval = NaN(n);

[R p] = chol(V); % right Cholesky factor
if p > 0 % fall back on "safe" method (doesn't ensure positive-definite result)
	R = inv(V);
	d = 1./sqrt(diag(R));
	R = -bsxfun(@times,bsxfun(@times,d,R),d');
else
	R = inv(R);
	R = bsxfun(@times,1./sqrt(sum(R.*conj(R),2)),R);
	R = -R*R';    % NOTE: so -R is positive-definite
end
R(1:n+1:n*n) = 1; % put 1s on the diagonal, so that each variable partially- correlates perfectly with itself

if nargout > 2 % calculate (approximate) p-values using Fisher z-transformation
	assert(nargin > 1 && ~isempty(N),'Need sample size for p-values');
	z = sqrt(N-n-1)*atanh(R);         % z-score for Fisher transformation (n-2 controlling variables)
	pval = 2*(1-normcdf(abs(z)));     % 2-tailed test
end
