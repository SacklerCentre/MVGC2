function [V,L,retries,iters] = cov_rand_alt2(n,r,corrmat,tol,maxretries,maxiters)

% Generate random positive-definite covariance matrix with specified generalised
% correlation. Sets up a multivariate normal random factor matrix, then adjusts
% the off-diagonal (using a binary chop) until within tolerance of specified
% generalised correlation.
%
% NOTE: generalised correlation is taken as r = sqrt(1-|R|^(2/n)) where R is the
% correlation matrix and n the number of dimensions. For n = 2 and small r, this
% is close to the Pearson correlation coefficient.
%
% WARNING: may well fail for gcorr ~ 2 or larger!
%
% n        - number of dimensions
% r        - generalised correlation: r = 0 yields zero correlation
% corrmat  - flag: return correlation matrix (default)
% tol      - numerical tolerance (default: sqrt(eps))
% maxiters - maximum iterations (default: 1000)
%
% V       - covariance matrix
% L       - Cholesky (left) factor: L*L' = V
% retries - number of retires required
% iters   - number of (binary chop) iterations required

assert(r >= 0,'generalised correlation must be non-negative');

gcorr = -log(1-r*r);

if nargin < 3 || isempty(corrmat),    corrmat    = true;      end
if nargin < 4 || isempty(tol),        tol        = sqrt(eps); end
if nargin < 5 || isempty(maxretries), maxretries = 1000;      end
if nargin < 5 || isempty(maxiters),   maxiters   = 1000;      end

odidx = 1:n*n; odidx(1:n+1:n*n) = []; % indices of off-diagonal entries

L = NaN;
retries = 0;
iters = 0;

if r < tol % diagonal cov matrix
    if corrmat
        L = eye(n);
        V = L;
    else
        L = diag(randn(n,1));
        V = L*L';
    end
    iters = 0;
    return
end

gotit = false;
fhi = 1;
maxhifacs  = 20;
for retries = 1:maxretries % retry until we are able to find a hi-factor large enough
    L = randn(n);
    F = L(odidx);
    for hfiters = 1:maxhifacs % find a hi-factor large enough
        L(odidx) = fhi*F;
        V = L*L';
		g = 2*(sum(log(diag(V)))-logdet(V))/n;
        if g > gcorr
            gotit = true;
            break
        end
        fhi = 2*fhi;
    end
    if gotit
        break
    end
end
if ~gotit % just not working - give up
    fprintf(2,'WARNING: ''cov_rand'' timed out on retries (is gcorr too large?)');
    V = NaN;
    return
end

gotit = false;
flo = 0;
for iters = 1:maxiters % binary chop
    f = (flo+fhi)/2;
    L(odidx) = f*F;
    V = L*L';
    g = 2*(sum(log(diag(V)))-logdet(V))/n;
    if     g > gcorr+tol
        fhi = f;
    elseif g < gcorr-tol
        flo = f;
    else
        gotit = true;
        break;
    end
end
if ~gotit % just not working - give up
    fprintf(2,'WARNING: ''cov_rand'' timed out on binary-chop (numerically unstable?)');
    V = NaN;
    return
end

L = chol(L*L','lower');
if corrmat
    L = diag(1./sqrt(sum(L.*L,2)))*L;
end
V = L*L';
