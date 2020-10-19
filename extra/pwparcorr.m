function [rho,pval] = pwparcorr(V,N)

% Calculate pairwise partial (Pearson) correlation coefficients from
% covariance matrix using inversion method. P-values are calculated
% for a 2-tailed test.
%
% V - sample covariance matrix (unbiased estimator)
% N - sample size (number of data observations)
%
% NOTE: The pairwise-conditional mutual information under Gaussian
% assumptions is -log(1-rho.*conj(rho)), with an asymptotic chi^2(1)
% distribution.

[n,n1] = size(V);
assert(ismatrix(V) && n1 == n,'Covariance matrix must be square');

rho  = NaN(n);
pval = NaN(n);

[V p] = chol(V);
assert(p == 0,'Covariance matrix not positive-definite');

V = inv(V);
rho = bsxfun(@times,1./sqrt(sum(V.*conj(V),2)),V);
rho = -rho*rho';     % NOTE: This means that -rho is positive-definite, but ...
rho(1:n+1:n*n) = 1;  % put 1s on the diagonal, so that each variable partially-
                     % correlates perfectly with itself (2I-rho is positive-definite)

if nargout > 1 % calculate (approximate) pval-values using Fisher z-transformation
	assert(nargin > 1 && ~isempty(N),'Need sample size for p-values');
	z = atanh(rho);                   % Fisher's z
	fac = sqrt(N-n-5);                % scale factor (there are n-2 controlling variables)
	pval = 2*(1-normcdf(fac*abs(z))); % 2-tailed test
end
