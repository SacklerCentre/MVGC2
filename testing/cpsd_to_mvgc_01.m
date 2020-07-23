function F = cpsd_to_mvgc(S,V,x,y,tol,maxiter)

if nargin < 5, tol     = []; end % cpsd_specfac default
if nargin < 6, maxiter = []; end % cpsd_specfac default

n = size(S,1);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);        % index of x in reduced quantities

F = NaN;

[~,VR,converged,relerr,niter] = cpsd_specfac(S(r,r,:),tol,maxiter)
if ~converged
    fprintf(2,'WARNING: spectral factorisation failed to converge in %d iterations (relative residual = %e)\n',niter,relerr);
end

F = logdet(VR(xr,xr)) - logdet(V(x,x));
