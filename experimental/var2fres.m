function [fres,frpow2,diff] = var2fres(A,V,tol,minpow2,maxpow2)

if nargin < 3 || isempty(tol), tol = 0; end

if abs(tol) < eps % i.e., zero

	% Calculate frequency resolution on the basis of autocorrelation decay to machine precision

	frpow2 = nextpow2(log(eps)/log(max(specnorm(A)))); % so that autocov decays to < eps
	fres = 2^frpow2;
	if nargout > 2 % want integral spectral check
		L = chol(V,'lower');
		LDV = 2*sum(log(diag(L)));
		H = var2trfun(A,fres);
		LDS = zeros(fres+1,1);
		for k = 1:fres+1 % over [0,pi]
			HLk = H(:,:,k)*L;
			LDS(k) = sum(log(diag(chol(HLk*HLk')))); % (log-determinant)/2
		end
		ILDS = sum(LDS(1:end-1)+LDS(2:end))/fres; % integrate frequency-domain logdet(S)
		diff = abs(ILDS-LDV);
	end

else

	% Calculate frequency resolution on the basis of spectral integral of logdet|S| = logdet|V|

	assert(nargin == 5,'Must supply minimum and maximum powers of 2');

	L = chol(V,'lower');
	LDV = 2*sum(log(diag(L)));
	failed = true;
	for frpow2 = minpow2:maxpow2
		fres = 2^frpow2;
		H = ss2trfun(A,fres);
		LDS = zeros(fres+1,1);
		for k = 1:fres+1 % over [0,pi]
			HLk = H(:,:,k)*L;
			LDS(k) = sum(log(diag(chol(HLk*HLk')))); % (log-determinant)/2
		end
		ILDS = sum(LDS(1:end-1)+LDS(2:end))/fres; % integrate frequency-domain logdet(S)
		diff = abs(ILDS-LDV);
		if diff <= tol
			failed = false;
			break;
		end
	end
	if failed, fres = NaN; end
end
