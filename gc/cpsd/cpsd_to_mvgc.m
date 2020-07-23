function F = cpsd_to_mvgc(S,x,y,tol,maxi,verb)

if nargin < 4, tol  = []; end % use cpsd_specfac default
if nargin < 5, maxi = []; end % use cpsd_specfac default
if nargin < 6, verb = 1;  end

n = size(S,1);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);        % index of x in reduced quantities

F = NaN;

[~,V,converged,relerr,niter] = cpsd_specfac(S,tol,maxi);
if converged
	if verb > 0
		fprintf('Full factorisation: iterations =%4d, rel. error = %e\n',niter,relerr);
	end
else
	if relerr > 0.5 % this has got to be a duff factorisation!
		fprintf(2,'ERROR: full spectral factorisation failed: iterations = %d, rel. error = %e\n',niter,relerr);
		return
	end
	fprintf(2,'WARNING: full spectral factorisation failed to converge: iterations = %d, rel. error = %e\n',niter,relerr);
end

[~,VR,converged,relerr,niter] = cpsd_specfac(S(r,r,:),tol,maxi);
if converged
	if verb > 0
		fprintf('Red. factorisation: iterations =%4d, rel. error = %e\n',niter,relerr);
	end
else
	if relerr > 0.5 % this has got to be a duff factorisation!
		fprintf(2,'ERROR: reduced spectral factorisation failed: iterations = %d, rel. error = %e\n',niter,relerr);
		return
	end
	fprintf(2,'WARNING: reduced spectral factorisation failed to converge: iterations = %d, rel. error = %e\n',niter,relerr);
end

F = logdet(VR(xr,xr)) - logdet(V(x,x));
