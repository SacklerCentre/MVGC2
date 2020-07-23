%-------------------------------------------------------------------------------
% Examine minimal bivariate processes
%-------------------------------------------------------------------------------

% Defaults

if ~exist('a',       'var'), a        = 0.8;     end % autoregressive parameter
if ~exist('c',       'var'), c        = 1.0;     end % causal feedback parameter
if ~exist('xi',      'var'), xi       = 1;       end % x residuals std. dev.
if ~exist('eta',     'var'), eta      = 1;       end % y residuals std. dev.
if ~exist('kappa',   'var'), kappa    = 0.5;     end % residuals correlation
if ~exist('m',       'var'), m        = 1000;    end % sequnce length
if ~exist('mtrans',  'var'), mtrans   = [];      end % transients to pre-truncate (empty for half sequnce length)
if ~exist('K',       'var'), K        = 20;      end % sequnce length
if ~exist('seed',    'var'), seed     = 0;       end % random seed
if ~exist('gpterm',  'var'), gpterm   = 'x-pdf'; end % Gnuplot temrinal

%-------------------------------------------------------------------------------

rstate = rng_seed(seed);
[x,y,epsx,epsy] = bvmin(cfb,a,c,xi,eta,kappa,m,mtrans);
rng_restore(rstate);

sx = std(x,1);
sy = std(y,1);

cxx = zeros(K+1,1);
cxy = zeros(K+1,1);
cyy = zeros(K+1,1);
for k = 0:K
	cxx(k+1) = mean(x(k+1:m).*x(1:m-k));
	cxy(k+1) = mean(x(k+1:m).*y(1:m-k));
	cyy(k+1) = mean(y(k+1:m).*y(1:m-k));
end

cx = cxx(1);
cy = cyy(1);

cxx = cxx/cx;
cxy = cxy/sqrt(cx*cy);
cyy = cyy/cy;

kk = (0:K)';

gp_qplot(kk,[cxx cyy cxy]);
