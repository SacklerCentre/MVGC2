function [fres,ierr] = var2fres(A,V,fastm,siparms)

if nargin < 2, V = []; end
if nargin < 3 || isempty(fastm),   fastm   = false;        end
if nargin < 4 || isempty(siparms), siparms = [1e-12,6,14]; end

assert(isvector(siparms) && length(siparms) == 3,'Spectral integration parameters must be a 3-vector [tolerance, fres pow2 min, fres pow2 max]');

tol     = siparms(1);
minpow2 = siparms(2);
maxpow2 = siparms(3);

if fastm

	% Calculate frequency resolution on the basis of autocorrelation decay

	frpow2 = nextpow2(log(eps)/log(specnorm(A))); % so that autocov decays to < eps
	if     frpow2 > maxpow2
		fprintf(2,'WARNING: Frequency resolution > 2^%d - defaulting to 2^%d\n',maxpow2,maxpow2);
		frpow2 = maxpow2;
	elseif frpow2 < minpow2
		fprintf(2,'WARNING: Frequency resolution < 2^%d - defaulting to 2^%d\n',minpow2,minpow2);
		frpow2 = minpow2;
	end
	fres = 2^frpow2;
	if nargout > 1 % want integral spectral check
		ierr = var_check_fres(A,V,fres);
	end

else

	% Calculate frequency resolution on the basis of spectral integral of logdet|S| = logdet|V|

	failed = true;
	for frpow2 = minpow2:maxpow2
		fres = 2^frpow2;
		ierr = var_check_fres(A,V,fres);
		if ierr <= tol
			failed = false;
			break;
		end
	end
	if failed
		fprintf(2,'WARNING: Frequency resolution > 2^%d - defaulting to 2^%d\n',maxpow2,maxpow2);
	end

end
