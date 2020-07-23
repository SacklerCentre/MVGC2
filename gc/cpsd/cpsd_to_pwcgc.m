function F = cpsd_to_pwcgc(S,tol,maxi,verb)

if nargin < 2, tol  = []; end % use cpsd_specfac default
if nargin < 3, maxi = []; end % use cpsd_specfac default
if nargin < 4, verb = 1;  end

n = size(S,1);

F = nan(n);

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
LDV = log(diag(V));

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

	[~,VR,converged,relerr,niter] = cpsd_specfac(S(r,r,:),tol,maxi);
	if converged
		if verb > 0
			fprintf('Red. factorisation, source %2d: iterations =%4d, rel. error = %e\n',y,niter,relerr);
		end
	else
		if relerr > 0.5 % this has got to be a duff factorisation!
			fprintf(2,'ERROR: reduced spectral factorisation failed, source %2d: iterations = %d, rel. error = %e\n',y,niter,relerr);
			continue
		end
		fprintf(2,'WARNING: reduced spectral factorisation failed to converge, source %2d: iterations = %d, rel. error = %e\n',y,niter,relerr);
	end

    F(r,y) = log(diag(VR))-LDV(r);
end
