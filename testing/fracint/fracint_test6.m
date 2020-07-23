tnet      = tnet5;  % connectivity network
%-------------------------------------------------------------------------------

pp        = 6;      % model order
rho       = 0.8;   % spectral radius
wvar      = 0.5;    % var coefficients decay weighting factor
rmi       = 0.5;    % residuals log-generalised correlation (multi-information)

fir       = 1000;
drange    = [0.1 0.5];

m         = 100000;
N         = 1;
mtrunc    = 10*m;

fs        = 200;
mtaper    = 1;
winfac    = [];

pregmode  = 'LWR';  % VAR model estimation regression mode ('OLS' or 'LWR')
psel      = 'LRT';  % model order selection ('ACT', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
pmax      = 40; % maximum model order for model order selection

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS' or 'LWR')

tstats    = 'dual'; % test statistic ('single', 'dual' or 'both')
alpha     = 0.05;   % significance level for Granger casuality significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fmin      = 0.1;
fmax      = fs/16;

q         = 10;

%-------------------------------------------------------------------------------

n = size(tnet,1);

d = drange(1)+(drange(2)-drange(1))*rand(n,1);

AA = var_rand(tnet,pp,rho,wvar);
VV = corr_rand(n,rmi);
infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

XX = varfima_to_tsdata(AA,[],{d fir true},VV,m,N,mtrunc);

h = zeros(n,1);
for i = 1:n
	h(i) = 0.5-Hurst_est(XX(i,:)',q);
end

% gp_qplot((1:m)',XX',[],[],[],[1 Inf]);

[SS,f,fres] = tsdata_to_cpsd(XX,mtaper,fs,winfac,[],-4,true,1); % autspectra only

fidx = find(f>=fmin&f<=fmax);
[a,b] = fitline(log(f(fidx)'),log(SS(:,fidx)));

fprintf('\nd = %s\n',num2str(d',' %6.4f'));
fprintf('h = %s\n',num2str(h',' %6.4f'));
fprintf('e = %s\n\n',num2str(-a'/2',' %6.4f'));
line1 = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead',fmin,fmin);
line2 = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead',fmax,fmax);
gpcmds = sprintf('set logs xy\n%s\n%s',line1,line2);
ff = [nan(fidx(1)-1,1);f(fidx);nan(fres-fidx(end)+1,1)];
SSf = exp(a*log(ff)'+b);
gp_qplot(f,[SS' SSf'],[],gpcmds,[],Inf);

return

X = gendiff(XX,{-d fir true});

[S,f,fres] = tsdata_to_cpsd(X,mtaper,fs,winfac,[],[],true,1); % autspectra only
gp_qplot(f,S',[],'set logs xy',[],Inf);


function [a,b] = fitline(x,y)

	n = size(x,2);
	xm  = sum(x,2)/n;
	ym  = sum(y,2)/n;
	xxm = sum(x.*x,2)/n;
	xym = sum(x.*y,2)/n;
	a   = (xym-xm.*ym)./(xxm-xm.*xm);
	b   = ym-xm.*a;

end
