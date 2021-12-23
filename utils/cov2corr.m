function [R,p,pval] = cov2corr(V,N)

% Calculate pairwise (Pearson) correlation matrix from covariance
% matrix. If the covariance matrix is an unbiased estimate from
% multivariate data, p-values are calculated for a 2-tailed test
% (i.e. for the null hypothesis of zero correlation); n(n-1)/2
% non-independent hypotheses must be assumed (see routine
% 'significance').
%
% The algorithm attempts to ensure a positive-definite result; if
% the covariance matrix is not positive-definite, we fall back on
% a "safe" method, for which the result may be inaccurate and not
% positive-definite; the flag p should thus always be checked.
%
% V - covariance matrix
% N - sample size (number of data observations)

n = size(V,1);
assert(issymmetric(V),'Covariance matrix must be square symmetric');

R  = NaN(n);
pval = NaN(n);

[R,p] = chol(V); % right Cholesky factor
if p > 0 % fall back on "safe" method (possibly inaccurate, and doesn't ensure positive-definite result)
	d = 1./sqrt(diag(V));
	R = symmetrise(bsxfun(@times,bsxfun(@times,d,V),d'));
else      % ensures positive-definite result
	R = bsxfun(@rdivide,R,sqrt(sum(R.*R)));
	R = R'*R;
end

if nargout > 2 % calculate (approximate) p-values using Fisher z-transformation
	assert(nargin > 1 && ~isempty(N),'Need sample size for p-values');
	RR = R; RR(1:n+1:n*n) = NaN;      % diagonal not relevant!
	z = sqrt(N-3)*atanh(RR);          % z-score for Fisher transformation
	pval = 2*(1-normcdf(abs(z)));     % 2-tailed test
end
