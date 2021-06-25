function fres = calc_fres(rho,tol)

if nargin < 2 || isempty(tol), tol = sqrt(eps); end

if rho < 1
	fres = 2^nextpow2(log(tol)/log(rho)); % so that autocorrelation decays to <= tol
else % assume number of observations
	fres = 2^nextpow2(rho);
end
