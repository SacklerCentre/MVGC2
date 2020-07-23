function [R,L,retries,iters] = corr_rand_alt(n,r,tol,maxretries,maxiters)

% Generate random correlation matrix with specified generalised
% correlation. Creates a uniform random orthogonal matrix and a random variance
% vector with independent chi^2(1) distribution, then adjusts the variance vector
% until within tolerance of specified generalised correlation.
%
% NOTE: generalised correlation is defined as r = |R| where R is the correlation
% matrix and n the number of dimensions. For n = 2, r = 1-rho^2 where rho is the
% Pearson correlation coefficient.
%
% n          - number of dimensions
% r          - generalised correlation: r = 1 yields zero correlation
% tol        - numerical tolerance (default: sqrt(eps))
% maxretries - maximum retries to find large enough correlation (default: 1000)
% maxiters   - maximum iterations for binary chop(default: 1000)
%
% R          - correlation matrix
% L          - Cholesky (left) factor: L*L' = R
% retries    - number of retries required
% iters      - number of binary chop iterations required

if nargin < 3 || isempty(tol),        tol        = sqrt(eps); end
if nargin < 4 || isempty(maxretries), maxretries = 1000;      end
if nargin < 5 || isempty(maxiters),   maxiters   = 1000;      end

assert(r >= 0 && r <= 1,'generalised correlation must be greater than 0 and less than or equal to 1');
if r < tol
	r = tol;
	fprintf(2,'WARNING: requested generalised correlation dangerously small: setting to %e\n',tol);
end

L       = NaN;
retries = 0;
iters   = 0;

if r > 1-tol
	R = eye(n);
	L = eye(n);
    return
end

% We calcluate a (positive-definite) covariance matrix with given generalised
% correlation, and finally convert it to a correlation matrix. For numerical
% efficiency, stability and precision, we work with the log form for generalised
% correlation, the multi-information:
%
% g = sum(log(diag(R)))-logdet(R))

gtarget = -log(r); % target value for g (see above)

% Find random orthogonal M and vector v of variances such that for R = M*diag(v)*M'
% (which will be pos-def), g is >= gtarget. The rationale is that adding the same
% constant c to all variances always decreases g, so we may use a binary chop to
% home in on gtarget. Note that since M is orthogonal, |R| = prod(v), so that we
% may calculate g efficiently as sum(log(diag(R))-log(v)).

gotit = false;
for retries = 0:maxretries
    [Q,R] = qr(randn(n));
    v = randn(n,1).^2;         % n independent chi^2(1)'s
    M = Q*diag(sign(diag(R))); % M orthogonal
    R = M*diag(v)*M';          % R is pos-def
    g = sum(log(diag(R))-log(v));
    if g >= gtarget
        gotit = true;
        break
    end
end
if ~gotit
    fprintf(2,'ERROR: ''corr_rand'' timed out on retries (is r too close to 0?)\n');
    R = NaN;
    returnR
end
D = diag(R);

% Got M and v (for the binary chop we just need the variances D of R and the
% unrotated variances v). Now set binary chop initial high value so that g <
% gtarget; start c at 1 and keep doubling

c = 1;
while g > gtarget
    g = sum(log(D+c) - log(v+c));
    c = 2*c;
end
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
    fprintf(2,'ERROR: ''corr_rand'' timed out on binary-chop (numerical instability?)\n');
    R = NaN;
    return
end

% We have a c that meets tolerance

R = M*diag(v+c)*M';

% Check R really is pos-def

[L,p] = chol(R,'lower');
if p ~= 0
    fprintf(2,'ERROR: ''corr_rand'' result not positive-definite\n');
    R = NaN;
    L = NaN;
    return
end

% Convert to correlation matrix

L = diag(1./sqrt(sum(L.*L,2)))*L;
R = L*L';
