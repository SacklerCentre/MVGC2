%-------------------------------------------------------------------------------
% Examine single-regression GC Type I errors
%-------------------------------------------------------------------------------

% Defaults

if ~exist('a',       'var'), a        = 0.6;     end % autoregressive parameter
if ~exist('c',       'var'), c        = 1.0;     end % causal feedback parameter
if ~exist('xi',      'var'), xi       = 1;       end % x residuals std. dev.
if ~exist('eta',     'var'), eta      = 1;       end % y residuals std. dev.
if ~exist('k',       'var'), k        = 0.7;     end % residuals correlation
if ~exist('m',       'var'), m        = 1000;    end % sequnce length
if ~exist('pmax',    'var'), pmax     = 10;      end % maximum AR model order
if ~exist('mtrans',  'var'), mtrans   = [];      end % transients to pre-truncate (empty for half sequnce length)
if ~exist('seed',    'var'), seed     = 0;       end % random seed
if ~exist('gpterm',  'var'), gpterm   = 'x-pdf'; end % Gnuplot temrinal

%-------------------------------------------------------------------------------

if isempty(mtrans), mtrans = round(m/2); end
mtot = m+mtrans;

rstate = rng_seed(seed);
u = randn(mtot,2);
rng_restore(rstate);
epsx = xi*u(:,1);
epsy = k*eta*u(:,1) + sqrt(1-k*k)*eta*u(:,2);

x = epsx;
y = epsy;
xy = x.*y;
for t = 2:mtot
	x(t) = x(t) + a*x(t-1) + c*xy(t-1);
end

%[mean(x)   (c*k*xi*eta)/(1-a)]
%[mean(xy)  k*xi*eta]

%gpcmds = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead front ls 2 lw 2',muthy,muthy);
%gp_qplot(t,x,[],gpcmds,gpterm);

%[moaic,mobic,mohqc,molrt] = tsdata_to_varmo([x y]',pmax,'LWR',[],true,gpterm);
%[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(x',    pmax,'LWR',[],true,gpterm);

X       = [x y]';
morder  = 1;
regmode = 'OLS';
tstats  = 'dual';
alpha   = 0.05;
mhtc    = 'FDR';

[A,V] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');

specnorm(A)

[F,stats] = var_to_pwcgc(A,V,tstats,X,regmode);
assert(~isbad(F,false),'GC estimation failed');

% Significance test (F- and likelihood ratio), adjusting for multiple hypotheses.

sigF  = significance(stats.(tstats).F.pval, alpha,mhtc);
sigLR = significance(stats.(tstats).LR.pval,alpha,mhtc);

display(F)
[stats.(tstats).F.pval stats.(tstats).LR.pval]
[sigF sigLR]

return

maxF = 1.1*nanmax(F(:));
pdata = {F,sigF,sigLR};
ptitle = {'PWCGC','sig(F)','sig(LR)'};
maxp = [maxF 1 1];
plot_gc(pdata,ptitle,[],maxp,gpterm);
