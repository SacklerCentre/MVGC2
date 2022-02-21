function [fres,ierr] = var2fres(A,V,siparms)

if nargin < 3 || isempty(siparms);

	% Calculate frequency resolution on the basis of autocorrelation decay

	hardmaxfrp2 = 16;
	frpow2 = nextpow2(log(eps)/log(specnorm(A))); % so that autocov decays to < eps
	if frpow2 > 2^hardmaxfrp2
		fprintf(2,'WARNING: Frequency > 2^%d - defaulting to 2^%d\n',hardmaxfrp2,hardmaxfrp2);
		frpow2 = 2^hardmaxfrp2;
	end
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
		ierr = abs(ILDS-LDV);
	end

else

	% Calculate frequency resolution on the basis of spectral integral of logdet|S| = logdet|V|

	assert(isvector(siparms) && length(siparms) == 3,'Spectral integration parameters must be a 3-vector [tolerance, fres pow2 min, fres pow2 max]');

	tol     = siparms(1);
	minpow2 = siparms(2);
	maxpow2 = siparms(3);

	L = chol(V,'lower');
	LDV = 2*sum(log(diag(L)));
	failed = true;
	for frpow2 = minpow2:maxpow2
		fres = 2^frpow2;
		H = var2trfun(A,fres);
		LDS = zeros(fres+1,1);
		for k = 1:fres+1 % over [0,pi]
			HLk = H(:,:,k)*L;
			LDS(k) = sum(log(diag(chol(HLk*HLk')))); % (log-determinant)/2
		end
		ILDS = sum(LDS(1:end-1)+LDS(2:end))/fres; % integrate frequency-domain logdet(S)
		ierr = abs(ILDS-LDV);
		if ierr <= tol
			failed = false;
			break;
		end
	end
	if failed, fres = NaN; end

end
