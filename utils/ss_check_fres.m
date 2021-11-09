function abserr = ss_check_fres(A,C,K,V,fres)

% Check that the log-determinant of the CPSD integrates to
% the log-determinant of the residuals covariance matrix at
% the given frequency resolution.

[L,cholp] = chol(V,'lower');
if cholp % show stopper
	abserr = NaN;
	return
end
LDV = 2*sum(log(diag(L)));

h     = fres+1;
[n,r] = size(C);
In    = eye(n);
Ir    = eye(r);
w     = exp(1i*pi*((0:fres)/fres));
LDS   = zeros(h,1);
for k = 1:h % over [0,pi]
	HLk = (In + C*((w(k)*Ir-A)\K))*L;
	LDS(k) = logdet(HLk*HLk');
end

abserr = abs(trapz(LDS)/fres - LDV);
