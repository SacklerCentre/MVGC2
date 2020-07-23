function [V,L] = cov_rand_alt1(n,icfac,corrmat)

% Generate random positive-definite covariance matrix
%
% n       - observation variable dimension
% icfac   - innovations correlation factor: a positive integer (as icfac gets
%           bigger, innovation correlations get smaller), or set to Inf for zero
%           correlation
% corrmat - flag: return correlation matrix (default)
%
% V       - covariance matrix
% L       - Cholesky (left) factor: L*L' = V
%
% V is generated as the sample covariance matrix of an n-dimensional
% uncorrelated Gaussian white noise sequence of length icfac*n. This ensures
% positive-definiteness. As ifac -> Inf, correlation between innovations -> 0.
% Setting icfac = Inf yields a diagonal V (no correlation). If the corrmat flag
% is set, V is normalised to variance 1; i.e. V is returned as a correlation
% matrix (the identity matrix in the case icfac = Inf).

assert(isscalar(icfac) && isnumeric(icfac) && icfac == floor(icfac) && icfac >= 1);

if nargin < 3 || isempty(corrmat), corrmat = true; end

if isinf(icfac) % diagonal covariance matrix
    if corrmat
        L = eye(n);
    else
        L = diag(randn(n,1));
    end
else
    L = randn(n,icfac*n)/sqrt(icfac*n); % sample covariance of n x n white noise of length icfac*n
    L = chol(L*L','lower');
    if corrmat
        L = diag(1./sqrt(sum(L.*L,2)))*L;
    end
end

V = L*L';
