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
if ~exist('p',       'var'), p        = 1;       end % AR model order
if ~exist('seed',    'var'), seed     = 0;       end % random seed
if ~exist('gpterm',  'var'), gpterm   = 'x-pdf'; end % Gnuplot temrinal

%-------------------------------------------------------------------------------

regmode = 'OLS';

F = zeros(S,2);

s = 1;
badm = 0;
badf = 0;
rstate = rng_seed(seed);
while true

	[x,y,epsx,epsy] = bvmin(cfb,a,c,xi,eta,kappa,m,mtrans);

	[A,V] = tsdata_to_var([x y]',p,regmode);
	if isbad(A) | specnorm(A) > 1-eps, badm = badm+1; continue; end

	% 2 -> 1
	[AR,VR] = tsdata_to_var(x',p,regmode);
	if isbad(AR) | abs(AR) > 1-eps, badf = badf+1; continue; end
	F(s,1) = logdet(VR)-logdet(V(1,1));

	% 1 -> 2
	[AR,VR] = tsdata_to_var(y',p,regmode);
	if isbad(AR) | abs(AR) > 1-eps, badf = badf+1; continue; end
	F(s,2) = logdet(VR)-logdet(V(2,2));

	s = s+1; if s > S, break; end
end
rng_restore(rstate);

pval = 1-chi2cdf((m-p)*F,p);

nbins = 20;
h = zeros(nbins,2);
h(1,1) = nnz(pval(:,1) <= 1/nbins);
h(1,2) = nnz(pval(:,2) <= 1/nbins);
for i = 2:nbins
	h(i,1) = nnz(pval(:,1) > (i-1)/nbins & pval(:,1) <= i/nbins);
	h(i,2) = nnz(pval(:,2) > (i-1)/nbins & pval(:,2) <= i/nbins);
end
cbin = (1:nbins)'/nbins - 1/(2*nbins);

gpstem = fullfile(tempdir,sprintf('%s_mode_%d_%d_m_%d_c_%4.2f',mfilename,cfb(1),cfb(2),m,c));
gp_write(gpstem,[cbin h/S]);
gp = gp_open(gpstem,gpterm,[2,0.6],14);
fprintf(gp,'datfile = "%s.dat"\n',gpstem);
%fprintf(gp,'set title "Type II error rate ($\\\\rho = %4.2f, F = %6.4f, M = %d, S = %d$)"\n',rho,Ft,M,S);
fprintf(gp,'set key top right Left rev\n');
fprintf(gp,'set yr [0:1]\n');
fprintf(gp,'set style fill solid 0.2\n');
fprintf(gp,'set xlabel "p-value"\n');
fprintf(gp,'set ylabel "frequency"\n');
fprintf(gp,'set multiplot layout 2,1\n');
fprintf(gp,'plot datfile using 1:2 w boxes t "$F(2 \\\\to 1)$"\n');
fprintf(gp,'plot datfile using 1:3 w boxes t "$F(1 \\\\to 2)$"\n');
fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm);
