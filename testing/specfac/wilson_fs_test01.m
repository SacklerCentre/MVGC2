%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('n',         'var'), n         = 7;        end % number of variables
if ~exist('p',         'var'), p         = 4;        end % AR model order
if ~exist('rho',       'var'), rho       = 0.8;      end % spectral norm
if ~exist('w',         'var'), w         = 1.0;      end % AR decay weighting parameter (empty for no weighting)
if ~exist('g',         'var'), g         = 1.0;      end % residuals multi-information (0 for zero residuals correlation, empty for uniform random)
if ~exist('fres',      'var'), fres      = 200;      end % frequency resolution
if ~exist('tol',       'var'), tol       = [];       end % tolerance for Wilson's algorithm
if ~exist('mseed',     'var'), mseed     = 0;        end % model random seed (0 to use current rng state)
if ~exist('gpterm',    'var'), gpterm    = 'x-pdf';  end % Gnuplot terminal type
if ~exist('verb',      'var'), verb      = 1;        end % verbosity level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate VAR model

rstate = rng_seed(mseed);
A = var_rand(n,p,rho,w);
V = corr_rand(n,g);
rng_restore(rstate);

% calculate spectrum, etc

[S,H] = var_to_cpsd(A,V,fres);

plot_cpsd(S);

% factorise

[H1,V1,ps,ps0,converged,relerr] = wilson_sf(S,1,tol);

h = fres+1;

for k = 1:h
	HL1k = H1(:,:,k)*chol(V,'lower');
	maxerrs(k) = maxabs(HL1k*HL1k'-S(:,:,k));
end
max(maxerrs)
