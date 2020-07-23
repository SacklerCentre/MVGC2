function f = cpsd_to_spwcgc(S,tol,maxi,verb)

if nargin < 2, tol  = []; end % use cpsd_specfac default
if nargin < 3, maxi = []; end % use cpsd_specfac default
if nargin < 4, verb = 1;  end

[n,~,h] = size(S);

[H,V,converged,relerr,niter] = cpsd_specfac(S,tol,maxi);
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

% for efficiency (in time if not memory) we pre-compute partial covariances

nr = n-1;
PL = zeros(nr,nr,n);
for x = 1:n
    w = [1:x-1 x+1:n]; % omit x
    PL(:,:,x) = chol(parcov(V,w,x),'lower'); % pre-compute the partial covariances for efficiency
end

BR = zeros(nr,nr,h);

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

	[HR,VR,converged,relerr,niter] = cpsd_specfac(S(r,r,:),tol,maxi);
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

	for k = 1:h
		BR(:,:,k) = inv(HR(:,:,k));
	end

    for xr = 1:n-1
        x  = r(xr);
        w = [1:x-1 x+1:n];  % omit x

        SR  = VR(xr,xr); % reduced model spectrum is flat!
        LSR = log(SR);
        PLx = PL(:,:,x);

        for k = 1:h
            HR = BR(xr,:,k)*H(r,w,k)*PLx;
            f(x,y,k) = LSR - log(SR-HR*HR');
        end
    end
end
