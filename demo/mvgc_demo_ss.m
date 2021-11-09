%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test data generation

ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial
fs        = 200;    % sample rate (Hz)

% Actual SS model generation parameters

nvars     = 5;      % number of variables
ssmoact   = 9;      % SS model order
rhoa      = 0.95;   % AR spectral radius
rmi       = 0.5;    % residuals log-generalised correlation (multi-information)
                    % g = -log|R|. g = 0 yields zero correlation,g = [] is uniform random
                    % on space of correlation matrices

% VAR model order estimation

varmosel  = 'AIC';  % VAR model order selection ('ACT', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
varmomax  = 32;     % maximum model order for VAR model order selection

% SS model order estimation

ssmosel   = 'SVC';  % SS model order selection ('ACT', 'SVC', or supplied numerical value)

% MVGC (frequency domain)

fres      = [];     % spectral MVGC frequency resolution (empty for automatic calculation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('seed',   'var'), seed     = 0;    end % random seed (0 for unseeded)
if ~exist('plotm',  'var'), plotm    = 0;    end % plot mode (figure number offset, or Gnuplot terminal string)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate random SS test data
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable X

% Seed random number generator.

rng_seed(seed);

% Generate random SS parameters in innovations form

[AA,CC,KK] = iss_rand(nvars,ssmoact,rhoa);

% Generate random residuals covariance (correlation) matrix.

VV = corr_rand(nvars,rmi);

% Report information on the generated SS model and check for errors.

infoo = ss_info(AA,CC,KK,VV);
assert(~infoo.error,'SS error(s) found - bailing out');

% Generate multi-trial SS time series data with normally distributed residuals
% for generated SS parameters and residuals covariance matrix.

ptic('*** ss_to_tsdata... ');
X = ss_to_tsdata(AA,CC,KK,VV,nobs,ntrials);
ptoc;

% Remove temporal mean and normalise by temporal variance.
% Not strictly necessary, but may help numerical stability
% if data has very large or very small values.

X = demean(X,true);

% VAR model order estimation (required for SS-SS estimation)

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.

ptic('\n*** tsdata_to_varmo... ');
if isnumeric(plotm), plotm = plotm+1; end
[varmoaic,varmobic,varmohqc,varmolrt] = tsdata_to_varmo(X,varmomax,'LWR',[],[],plotm);
ptoc;

% Select and report VAR model order.

varmo = moselect(sprintf('VAR model order selection (max = %d)',varmomax),varmosel,'AIC',varmoaic,'BIC',varmobic,'HQC',varmohqc,'LRT',varmolrt);
assert(varmo > 0,'selected zero model order! GCs will all be zero!');
if varmo >= varmomax, fprintf(2,'*** WARNING: selected VAR maximum model order (may have been set too low)\n'); end

% SS model order estimation

pf = 2*varmo; % Bauer recommends 2 x VAR model order

ptic('\n*** tsdata_to_ssmo... ');
if isnumeric(plotm), plotm = plotm+1; end
[ssmosvc,ssmomax] = tsdata_to_ssmo(X,pf,plotm);
ptoc;

% Select and report SS model order.

ssmo = moselect(sprintf('SS model order selection (max = %d)',ssmomax),ssmosel,'ACT',ssmoact,'SVC',ssmosvc);

assert(ssmo > 0,'selected zero model order! GCs will all be zero!');
if ssmo >= ssmomax, fprintf(2,'*** WARNING: selected SS maximum model order (may have been set too low)\n'); end

% SS model estimation

% Estimate SS model order and model paramaters using CCA SS-SS algorithm

[A,C,K,V] = tsdata_to_ss(X,pf,ssmo);

% Report information on the estimated SS, and check for errors.

info = ss_info(A,C,K,V);
assert(~info.error,'SS error(s) found - bailing out');

% Granger causality calculation: time domain

% Estimated time-domain pairwise-conditional Granger causalities

ptic('*** ss_to_pwcgc... ');
F = ss_to_pwcgc(A,C,K,V);
ptoc;
assert(~isbad(F,false),'GC estimation failed');

% NOTE: we don't have an analytic (asymptotic) distribution for the statistic, so no significance testing here!

% For comparison, we also calculate the actual pairwise-conditional causalities

ptic('*** ss_to_pwcgc... ');
FF = ss_to_pwcgc(AA,CC,KK,VV);
ptoc;
assert(~isbad(FF,false),'GC calculation failed');

% Plot time-domain causal graph

maxF = 1.1*max(nanmax(F(:),nanmax(FF(:))));
if isnumeric(plotm), plotm = plotm+1; end
plot_gc({FF,F},{'PWCGC (actual)','PWCGC (estimated)'},[],[maxF maxF],plotm);

%% Granger causality estimation: frequency domain

% Calculate spectral pairwise-conditional causalities resolution from SS model
% parameters. If not specified, we set the frequency resolution to something
% sensible. Warn if resolution is very large, as this may cause problems.

if isempty(fres)
	maxfres = 2^14; % adjust to taste
	fres = max(info.fres,infoo.fres);
	if fres > maxfres
		fprintf(2,'\nWARNING: esitmated frequency resolution %d exceeds maximum; setting to %d' ,fres,maxfres);
		fres = maxfres;
	else
		fprintf('\nUsing frequency resolution %d',fres);
	end
end
fabserr = ss_check_fres(A,C,K,V,fres);
fprintf(' (absolute integration error = %e)\n',fabserr);

ptic(sprintf('\n*** ss_to_spwcgc (at frequency resolution = %d)... ',fres));
f = ss_to_spwcgc(A,C,K,V,fres);
ptoc;
assert(~isbad(f,false),'spectral GC estimation failed');

% For comparison, we also calculate the actual pairwise-conditional spectral causalities

ptic(sprintf('*** ss_to_spwcgc (at frequency resolution = %d)... ',fres));
ff = ss_to_spwcgc(AA,CC,KK,VV,fres);
ptoc;
assert(~isbad(ff,false),'spectral GC calculation failed');

% Get frequency vector according to the sampling rate.

freqs = sfreqs(fres,fs);

% Plot spectral causal graphs.

if isnumeric(plotm), plotm = plotm+1; end
plot_sgc({ff,f},freqs,'Spectral Granger causalities (blue = actual, red = estimated)',[],plotm);

% Granger causality calculation: frequency domain -> time-domain

% Check that spectral causalities average (integrate) to time-domain
% causalities. Note that this may occasionally fail if a certain condition
% on the VAR parameters is not satisfied (Geweke 1982).

Fint = bandlimit(f,3); % integrate spectral MVGCs (frequency is dimension 3 of CPSD array)

fprintf('\n*** GC spectral integral check... ');
rr = abs(F-Fint)./(1+abs(F)+abs(Fint)); % relative residuals
mrr = max(rr(:));                       % maximum relative residual
if mrr < 1e-5
    fprintf('PASS: max relative residual = %.2e\n',mrr);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',mrr);
end
