tnet      = tnet9;  % connectivity network
pp        = 6;      % model order
rho       = 0.8;    % spectral radius
wvar      = 0.5;    % var coefficients decay weighting factor
rmi       = 0.5;    % residuals log-generalised correlation (multi-information)

fir       = 1000;
drange    = [0.4 0.5];

m         = 1000;
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

%-------------------------------------------------------------------------------

n = size(tnet,1);

d = drange(1)+(drange(2)-drange(1))*rand(n,1);
fprintf('\nd = %s\n',num2str(d',' %6.4f'));
fiparms = {d fir true};

AA = var_rand(tnet,pp,rho,wvar);
VV = corr_rand(n,rmi);
infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

XX = varfima_to_tsdata(AA,[],fiparms,VV,m,N,mtrunc);

gp_qplot((1:m)',XX',[],[],[],[1 Inf]);

[S,f,fres] = tsdata_to_cpsd(XX,mtaper,fs,winfac,[],[],true,1); % autspectra only
gp_qplot(f,S',[],'set logs xy',[],Inf);

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo_gp(XX,pmax,pregmode);

p = moselect(psel,'actual',pp,'AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
assert(p > 0,'selected zero model order! GCs will all be zero!');
if p >= pmax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end

[A,V] = tsdata_to_var(XX,p,regmode);
info = var_info(A,V);
assert(~info.error,'VAR error(s) found - bailing out');

[F,stats] = var_to_pwcgc(A,V,tstats,XX,regmode);
assert(~isbad(F,false),'GC estimation failed');

sigF  = significance(stats.(tstats).F.pval, alpha,mhtc);
sigLR = significance(stats.(tstats).LR.pval,alpha,mhtc);

FF = var_to_pwcgc(AA,VV);
assert(~isbad(FF,false),'GC calculation failed');

disp('actual GCs =');      disp(FF)
disp('estimated GCs = ='); disp(F)

tnet(1:n+1:n*n) = NaN;
sigact = tnet( ~isnan(tnet ));
sigf   = sigF( ~isnan(sigF ));
siglr  = sigLR(~isnan(sigLR));

fprintf('F  errors = %d\n',nnz(abs(sigact-sigf)  > eps));
fprintf('LR errors = %d\n',nnz(abs(sigact-siglr) > eps));
