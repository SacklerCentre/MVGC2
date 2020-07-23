function [V,L,retries,iters] = cov_rand_alt3(n,r,tol,maxretries,maxiters)

% Generate random correlation matrix with specified generalised
% correlation. Creates a uniform random orthogonal matrix and a random variance
% vector with independent chi^2(1) distribution, then adjusts the variance vector
% until within tolerance of specified generalised correlation.
%
% NOTE: generalised correlation is defined as r = |R| where R is the correlation
% matrix and n the number of dimensions. For n = 2 and small r, r = 1-rho^2 where
% rho is the Pearson correlation coefficient.
%
% WARNING: may fail for r close to 0!
%
% n          - number of dimensions
% r          - generalised correlation: r = 0 yields zero correlation
% tol        - numerical tolerance (default: sqrt(eps))
% maxretries - maximum retries to find large enough correlation (default: 1000)
% maxiters   - maximum iterations for binary chop(default: 1000)
%
% V          - correlation matrix
% L          - Cholesky (left) factor: L*L' = V
% retries    - number of retires required
% iters      - number of (binary chop) iterations required

assert(r >= 0 && r <= 1,'generalised correlation must lie between 0 and 1');

if nargin < 3 || isempty(tol),        tol        = sqrt(eps); end
if nargin < 4 || isempty(maxretries), maxretries = 1000;      end
if nargin < 5 || isempty(maxiters),   maxiters   = 1000;      end

L       = NaN;
retries = 0;
iters   = 0;

if r > 1-tol % diagonal cov matrix
	L = eye(n);
	V = eye(n);
    return
end

% We calcluate a (positive-definite) covariance matrix with given generalised
% correlation. For numerical efficiency, stability and precision, we work with the
% log form for generalised correlation (multi-information):
%
% g = sum(log(diag(V)))-logdet(V))

gtarget = -log(r) % target value for g (see above)

% Find random orthogonal M and vector v of variances such that for V =
% M*diag(v)*M' (which will be pos-def), g is >= gtarget. The rationale is that
% adding the same constant c to all variances always decreases g, so we may use
% a binary chop to home in on gtarget.
%
% Note that
%
% g = sum(log(diag(V))-log(v))
%
% is efficient to calculate.

gotit = false;
for retries = 0:maxretries
    [Q,R] = qr(randn(n));
    v = randn(n,1).^2;         % n independent chi^2(1)'s
    M = Q*diag(sign(diag(R))); % M orthogonal
    V = M*diag(v)*M';          % V is pos-def
    g = sum(log(diag(V))-log(v))
    if g >= gtarget
        gotit = true;
        break
    end
end
if ~gotit
    fprintf(2,'ERROR: ''cov_rand'' timed out on retries (is r too close to 1?)\n');
    V = NaN;
    return
end
D = diag(V);

% Got M and v (for the binary chop we just need the variances D of V and the
% unrotated variances v). Now set binary chop initial high value so that g <
% gtarget; start c at 1 and keep doubling

c = 1;
while g > gtarget
    g = sum(log(D+c) - log(v+c));
    c = 2*c;
end
g
c
chi = c;
clo = 0;
ghi = gtarget+tol;
glo = gtarget-tol;

% Do binary chop

gotit = false;
for iters = 1:maxiters % binary chop
    c = (clo+chi)/2;
    g = sum(log(D+c) - log(v+c));
    if     g < glo % g too small, set chi to current
        chi = c;
    elseif g > ghi % g too big, set clo to current
        clo = c;
    else
        gotit = true;
        break;
    end
end
if ~gotit % this shouldn't really happen :-/
    fprintf(2,'ERROR: ''cov_rand'' timed out on binary-chop (numerical instability?)\n');
    V = NaN;
    return
end

% We have a c that meets tolerance

V = M*diag(v+c)*M';

% Check V really is pos-def

[L,p] = chol(V,'lower');
if p ~= 0
    fprintf(2,'ERROR: ''cov_rand'' result not positive-definite\n');
    V = NaN;
    L = NaN;
    return
end

% Convert to correlation matrix

L = diag(1./sqrt(sum(L.*L,2)))*L;
V = L*L';
