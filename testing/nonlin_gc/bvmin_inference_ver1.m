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
if ~exist('mtrans',  'var'), mtrans   = [];      end % transients to pre-truncate (empty for half sequnce length)
if ~exist('p',       'var'), p        = 1;       end % AR model order
if ~exist('seed',    'var'), seed     = 0;       end % random seed
if ~exist('gpterm',  'var'), gpterm   = 'x-pdf'; end % Gnuplot temrinal

%-------------------------------------------------------------------------------

regmode = 'OLS';
tstats  = 'dual';
alpha   = 0.05;
mhtc    = 'FDR';

FX = zeros(S,2);
%FF = zeros(S,2);

s = 1;
badm = 0;
badf = 0;
rstate = rng_seed(seed);
while true

	[x,y,epsx,epsy] = bvmin(cfb,a,c,xi,eta,kappa,m,mtrans);

	[A,V] = tsdata_to_var([x y]',p,regmode);
	if isbad(A) || specnorm(A) > 1-eps, badm = badm+1; continue; end

	% 2 -> 1
	[AR,VR] = tsdata_to_var(x',p,regmode);
	if isbad(AR) || abs(AR) > 1-eps, badf = badf+1; continue; end
	FX(s,1) = logdet(VR)-logdet(V(1,1));
	%FF(s,1) = VR/V(1,1)-1;

	% 1 -> 2
	[AR,VR] = tsdata_to_var(y',p,regmode);
	if isbad(AR) || abs(AR) > 1-eps, badf = badf+1; continue; end
	FX(s,2) = logdet(VR)-logdet(V(2,2));
	%FF(s,2) = VR/V(2,2)-1;

%{
	% 2 -> 1
	[~,v1,rep] = var2riss(A,V,2,1);
	if sserror(rep), badf = badf+1; continue; end
	F(s,1) = v1/V(1,1)-1;

	% 1 -> 2
	[~,v2,rep] = var2riss(A,V,1,2);
	if sserror(rep), badf = badf+1; continue; end
	F(s,2) = v2/V(2,2)-1;
%}

	s = s+1; if s > S, break; end
end
rng_restore(rstate);

d  = p;       % df
pvalx = 1-chi2cdf((m-p)*FX,d);

%d2 = m-3*p-1; % F df2
%K  = d2/d;    % F scaling factor
%pvalf = 1-fcdf(K*FF,d,d2);

nbins = 20;
h = zeros(nbins,2);
h(1,1) = nnz(pvalx(:,1) <= 1/nbins);
h(1,2) = nnz(pvalx(:,2) <= 1/nbins);
for i = 2:nbins
	h(i,1) = nnz(pvalx(:,1) > (i-1)/nbins & pvalx(:,1) <= i/nbins);
	h(i,2) = nnz(pvalx(:,2) > (i-1)/nbins & pvalx(:,2) <= i/nbins);
end
cbin = (1:nbins)'/nbins - 1/(2*nbins);

gpstem = fullfile(tempdir,sprintf('%s_%d',mfilename,cfb));
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
