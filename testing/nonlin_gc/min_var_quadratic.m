%-------------------------------------------------------------------------------
% Examine single-regression GC Type I errors
%-------------------------------------------------------------------------------

% Defaults

if ~exist('a',       'var'), a        = 0.5;     end % autoregressive parameter
if ~exist('c',       'var'), c        = 1.0;     end % causal feedback parameter parameter
if ~exist('xi',      'var'), xi       = 1;       end % x std. dev.
if ~exist('eta',     'var'), eta      = 1;       end % y std. dev.
if ~exist('m',       'var'), m        = 1000;    end % sequnce length
if ~exist('pmax',    'var'), pmax     = 40;      end % maximum AR model order
if ~exist('mtrans',  'var'), mtrans   = [];      end % transients to pre-truncate (empty for half sequnce length)
if ~exist('seed',    'var'), seed     = 0;       end % random seed
if ~exist('gpterm',  'var'), gpterm   = 'x-pdf'; end % Gnuplot temrinal

%-------------------------------------------------------------------------------

if isempty(mtrans), mtrans = round(m/2); end
mtot = m+mtrans;

rstate = rng_seed(seed);
epsx = xi*randn(mtot,1);
epsy = eta*randn(mtot,1);
rng_restore(rstate);

x = epsx;
y = epsy;
y2 = y.*y;
for t = 2:mtot
	x(t) = x(t) + a*x(t-1) + c*y2(t-1);
end

xL   = x (mtrans:mtot-1);
yL   = y (mtrans:mtot-1);
y2L  = y2(mtrans:mtot-1);

epsx = epsx(mtrans+1:mtot);
epsy = epsy(mtrans+1:mtot);
x    = x   (mtrans+1:mtot);
y    = y   (mtrans+1:mtot);
y2   = y2  (mtrans+1:mtot);
t = (1:m)';

a2 = a*a;
c2 = c*c;
xi2 = xi*xi;
eta2 = eta*eta;
eta4 = eta2*eta2;

[mean(xL.*y2) c*eta4/(1-a)]
[mean(x.*x) (xi2 + ((3-a)/(1-a))*c2*eta4)/(1-a2)]

muemp = mean(x);
muthy = (c/(1-a))*eta2;
[muemp muthy]

sig2emp = var(x);
sig2thy = (xi2 + 2*c2*eta4)/(1-a2);
[sig2emp sig2thy]

gpcmds = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead ls 2',muthy,muthy);
%gp_qplot(t,x,[],gpcmds,gpterm);

%[moaic,mobic,mohqc,molrt] = tsdata_to_varmo([x y]',pmax,'LWR',[],true,gpterm);

X       = [x y]';
morder  = 1;
regmode = 'LWR';
tstats  = 'dual';
alpha   = 0.05;
mhtc    = 'FDR';

[A,V] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');

[F,stats] = var_to_pwcgc(A,V,tstats,X,regmode);
assert(~isbad(F,false),'GC estimation failed');

% Significance test (F- and likelihood ratio), adjusting for multiple hypotheses.

sigF  = significance(stats.(tstats).F.pval, alpha,mhtc);
sigLR = significance(stats.(tstats).LR.pval,alpha,mhtc);

display(F)
pvalF = stats.(tstats).F.pval
sigF
pvalLR = stats.(tstats).LR.pval
sigLR

maxF = 1.1*nanmax(F(:));
pdata = {F,sigF,sigLR};
ptitle = {'PWCGC','sig(F)','sig(LR)'};
maxp = [maxF 1 1];
plot_gc(pdata,ptitle,[],maxp,gpterm);
