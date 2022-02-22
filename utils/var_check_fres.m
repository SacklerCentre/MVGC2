function ierr = var_check_fres(A,V,fres)

% Check that the log-determinant of the CPSD integrates to
% the log-determinant of the residuals covariance matrix at
% the given frequency resolution.

[L,cholp] = chol(V,'lower');
if cholp % show stopper
	ierr = NaN;
	return
end
LDV = 2*sum(log(diag(L)));
H = var2trfun(A,fres);
LDS = zeros(fres+1,1);
for k = 1:fres+1 % over [0,pi]
	HLk = H(:,:,k)*L;
	LDS(k) = sum(log(diag(chol(HLk*HLk')))); % (log-determinant)/2
end
ILDS = sum(LDS(1:end-1)+LDS(2:end))/fres; % integrate frequency-domain logdet(S)
ierr = abs(ILDS-LDV);
