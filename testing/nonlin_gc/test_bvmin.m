%-------------------------------------------------------------------------------
% Examine minimal bivariate processes
%-------------------------------------------------------------------------------

% Defaults

if ~exist('a',       'var'), a        = 0.8;     end % autoregressive parameter
if ~exist('c',       'var'), c        = 1.0;     end % causal feedback parameter
if ~exist('xi',      'var'), xi       = 2;       end % x residuals std. dev.
if ~exist('eta',     'var'), eta      = 5;       end % y residuals std. dev.
if ~exist('kappa',   'var'), kappa    = 0.4;     end % residuals correlation
if ~exist('m',       'var'), m        = 10000;   end % sequnce length
if ~exist('pmax',    'var'), pmax     = 20;      end % maximum AR model order
if ~exist('mtrans',  'var'), mtrans   = [];      end % transients to pre-truncate (empty for half sequnce length)
if ~exist('seed',    'var'), seed     = 0;       end % random seed
if ~exist('gpterm',  'var'), gpterm   = 'x-pdf'; end % Gnuplot temrinal

%-------------------------------------------------------------------------------

rstate = rng_seed(seed);
[x,y,epsx,epsy] = bvmin(cfb,a,c,xi,eta,kappa,m,mtrans);
rng_restore(rstate)

t = (1:m)';

gp_qplot(t,x);

X = [x y]';

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo([x y]',pmax,'LWR',[],[],'x11');
