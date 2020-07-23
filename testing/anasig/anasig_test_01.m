
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('nobs',  'var'), nobs  = 1001;    end
if ~exist('varmo', 'var'), varmo = 6;       end
if ~exist('mu',    'var'), mu    = 4.0;     end
if ~exist('rho',   'var'), rho   = 0.9;     end
if ~exist('wvar',  'var'), wvar  = 0.5;     end
if ~exist('rmi',   'var'), rmi   = 0.5;     end
if ~exist('fs',    'var'), fs    = 250;     end
if ~exist('fres',  'var'), fres  = [];      end
if ~exist('mosel', 'var'), mosel = 'HQC';   end
if ~exist('momax', 'var'), momax = 3*varmo; end
if ~exist('alpha', 'var'), alpha = 0.05;    end
if ~exist('mhtc',  'var'), mhtc  = 'FDRD';  end
if ~exist('seed',  'var'), seed  = 0;       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed random number generator.

rng_seed(seed);

% Generate random VAR coefficients for test network.

AA = var_rand(tnet5,varmo,rho,wvar);
nvars = size(AA,1); % number of variables

% Generate random residuals covariance (in fact correlation) matrix.

VV = corr_rand(nvars,rmi);
infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

% Generate VAR time series data

X = varfima_to_tsdata(AA,[],[],VV,nobs) + mu;

% Calculate analytic signal

Y = abs(hilbert(X')');

% Time stamp (seconds)

t = (0:nobs-1)'/fs;

% Remove temporal means

X = demean(X);
Y = demean(Y);

% Model order estimation

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(Y,momax,'LWR',[],[],1);

% Select and report VAR model order.

morder = moselect(sprintf('VAR model order selection (max = %d)',momax), mosel,'ACT',varmo,'AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
assert(morder > 0,'selected zero model order! GCs will all be zero!');
if morder >= momax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end

% Estimate VAR model of selected order from data.

[A,V] = tsdata_to_var(Y,morder,'LWR');

assert(~isbad(A),'VAR estimation failed - bailing out');

% Report information on the estimated VAR, and check for errors.

info = var_info(A,V);
assert(~info.error,'VAR error(s) found - bailing out');

% Granger causality calculation: time domain

[F,pval] = var_to_pwcgc(A,V,Y,'LWR');
assert(~isbad(F,false),'GC estimation failed');

% Significance test (F-test and likelihood ratio), adjusting for multiple hypotheses.

sigFT = significance(pval.FT,alpha,mhtc);
sigLR = significance(pval.LR,alpha,mhtc);

% For comparison, we also calculate the actual pairwise-conditional causalities

FF = var_to_pwcgc(AA,VV);
assert(~isbad(FF,false),'GC calculation failed');

% Plot time-domain causal graph and significance.

maxF = 1.1*max(nanmax(F(:),nanmax(FF(:))));
pdata = {FF,F;sigFT,sigLR};
ptitle = {'PWCGC (actual)','PWCGC (estimated)'; 'F-test','LR test'};
maxp = [maxF maxF;1 1];
plot_gc(pdata,ptitle,[],maxp,2);

% Relative error

relerr = (FF-F)/maxabs(FF);
fprintf('\nRelative error =\n');
disp(relerr);

%% Granger causality estimation: frequency domain

if isempty(fres)
    fres = 2^nextpow2(max(info.acdec,infoo.acdec)); % alternatively, fres = 2^nextpow2(nobs);
	fprintf('\nUsing frequency resolution %d\n',fres);
end
if fres > 10000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

f = var_to_spwcgc(A,V,fres);
assert(~isbad(f,false),'spectral GC estimation failed');

% For comparison, we also calculate the actual pairwise-conditional spectral causalities

ff = var_to_spwcgc(AA,VV,fres);
assert(~isbad(ff,false),'spectral GC calculation failed');

% Get frequency vector according to the sampling rate.

freqs = sfreqs(fres,fs);

% Plot spectral causal graphs.

plot_sgc({ff,f},freqs,'Spectral Granger causalities (blue = actual, red = estimated)',3);

% Granger causality calculation: frequency domain -> time-domain

Fint = bandlimit(f,3); % integrate spectral MVGCs (frequency is dimension 3 of CPSD array)

fprintf('\n*** GC spectral integral check... ');
rr = abs(F-Fint)./(1+abs(F)+abs(Fint)); % relative residuals
mrr = max(rr(:));                       % maximum relative residual
if mrr < 1e-6
    fprintf('PASS: max relative residual = %.2e\n',mrr);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',mrr);
end
