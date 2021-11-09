function abserr = var_check_fres(A,V,fres)

% Check that the log-determinant of the CPSD integrates to
% the log-determinant of the residuals covariance matrix at
% the given frequency resolution.

[L,cholp] = chol(V,'lower');
if cholp % show stopper
	abserr = NaN;
	return
end
LDV = 2*sum(log(diag(L)));

h   = fres+1;
n = size(A,1);
LDS = zeros(h,1);
AFT = fft(cat(3,eye(n),-A),2*fres,3); % over [0,2*pi)
for k = 1:h
	HLk = AFT(:,:,k)\L;
	LDS(k) = logdet(HLk*HLk');
end

abserr = abs(trapz(LDS)/fres - LDV);
