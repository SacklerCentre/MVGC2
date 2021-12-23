function [R,p,pval] = cov2pcorr(V,N)

% Calculate pairwise partial (Pearson) correlation matrix from
% covariance matrix using inversion method. If the covariance
% matrix is an unbiased estimate from multivariate data, p-values
% are calculated for a 2-tailed test (i.e. for the null hypothesis
% of zero partial correlation); n(n-1)/2 non-independent hypotheses
% must be assumed (see routine 'significance').
%
% The algorithm attempts to ensure a positive-definite result; if
% the covariance matrix is not positive-definite, we fall back on
% a "safe" method, for which the result may be inaccurate and not
% positive-definite; the flag p should thus always be checked.
%
% V - covariance matrix
% N - sample size (number of data observations)
%
% NOTE 1: For partial correlation calculated from multivariate data,
% this routine may be less accurate than Matlab's 'partialcorr'
% function - but it is definitely more efficient!
%
% NOTE 2: The pairwise-conditional mutual information under Gaussian
% assumptions is -log(1-R.*conj(R)), with an asymptotic chi^2(1)
% distribution.

n = size(V,1);
assert(issymmetric(V),'Covariance matrix must be square symmetric');

R  = NaN(n);
pval = NaN(n);

[R p] = chol(V); % right Cholesky factor
if p > 0 % fall back on "safe" method (possibly inaccurate, and doesn't ensure positive-definite result)
	warning('off','MATLAB:nearlySingularMatrix'); % if not positive-definite, expect problems anyway!
	R = inv(V);
	warning('on','MATLAB:nearlySingularMatrix');
	d = 1./sqrt(diag(R));
	R = symmetrise(-bsxfun(@times,bsxfun(@times,d,R),d'));
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
