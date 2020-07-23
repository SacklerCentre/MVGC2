function f = cpsd_to_smvgc(S,x,y,tol,maxi,verb)

if nargin < 4, tol  = []; end % use cpsd_specfac default
if nargin < 5, maxi = []; end % use cpsd_specfac default
if nargin < 6, verb = 1;  end

[n,~,h] = size(S);
fres = h-1;

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
w = [y z];
xr = 1:length(x);

f = nan(1,h);

[H,V,converged,relerr,niter,~,L] = cpsd_specfac(S,tol,maxi);
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
PL = chol(parcov(V,w,x),'lower');

% Note: in theory we shouldn't have to take the real part of the determinants of
% the (Hermitian, positive-definite) matrices in the calculation of the f(k),
% since they should be real. However, Matlab's det() will sometimes return a
% tiny imaginary part.

if isempty(z) % unconditional (note: does not require reduced model)

    for k = 1:h
        HL   = H(x,:,k)*L;
        SR   = HL*HL';
        HR   = H(x,y,k)*PL;
        f(k) = logdet(SR) - logdet(SR-HR*HR');
    end

else % conditional

	[HR,VR,converged,relerr,niter] = cpsd_specfac(S(r,r,:),tol,maxi);
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

	nr = length(r);
	BR = zeros(nr,nr,h);
	for k = 1:h
		BR(:,:,k) = inv(HR(:,:,k));
	end

    SR    = VR(xr,xr); % reduced model spectrum is flat!
    LDSR  = logdet(SR);
    for k = 1:h
        HR   = BR(xr,:,k)*H(r,w,k)*PL;
        f(k) = LDSR - logdet(SR-HR*HR');
    end

end
