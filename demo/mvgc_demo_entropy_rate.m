%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test data generation parameters[NOTE: omit if supplying your own time series data]

defvar('nchans',    5     ); % number of channels
defvar('nepochs',   10    ); % number of epochs
defvar('nobs',      1000  ); % number of observations per epoch
defvar('ssmoact',   9     ); % ISS model order
defvar('rho',       0.95  ); % AR spectral radius
defvar('rmi',       1     ); % residuals log-generalised correlation (multi-information)
defvar('seed',      0     ); % random seed (0 for unseeded)

% Other parameters

defvar('fs',        200   ); % sample rate (Hz)
defvar('ernorm',    false ); % "normalise" entropy rate (i.e., subtract process entropy)
defvar('varmomax',  32    ); % maximum model order for VAR model order selection (required for ISS estimation)
defvar('fres',      []    ); % frequency resolution  for spectral entropy rate (empty for automatic calculation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random ISS test data
%
% NOTE: This is where you would read in your own time series data; it should
% be assigned to the variable X as [channels x observations x epochs]; then set
%
%     [nchans,nobs,nepochs] = size(X);
%
% and skip the following code down to the next section.

% Seed random number generator

rng_seed(seed);

% Generate random ISS parameters

[AA,CC,KK] = iss_rand(nchans,ssmoact,rho);

% Generate random residuals covariance (correlation) matrix

VV = corr_rand(nchans,rmi);

% Report information on the generated ISS model

infoo = ss_info(AA,CC,KK,VV);
assert(~infoo.error,'ISS error(s) found - bailing out');

% Generate multi-epoch ISS time series data with normally distributed residuals

ptic('*** ss_to_tsdata... ');
X = ss_to_tsdata(AA,CC,KK,VV,nobs,nepochs);
ptoc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VAR model order estimation (required for ISS estimation)

ptic('\n*** tsdata_to_varmo... ');
varmo = tsdata_to_varmo(X,varmomax,'LWR',[],false,1);
ptoc;
assert(varmo > 0,'selected zero VAR model order (non-stationary data or other badness?)');
if varmo >= varmomax, fprintf(2,'*** WARNING: selected maximum VAR model order (''varmomax'' set too low?)\n'); end
fprintf('\nVAR model order (AIC) = %d\n',varmo);

% ISS model order estimation

pf = 2*varmo; % past/future horizon (Bauer recommends 2 x VAR AIC model order)

ptic('\n*** tsdata_to_ssmo... ');
[ssmo,ssmomax] = tsdata_to_ssmo(X,pf,2);
ptoc;

% Check and report ISS model order

assert(ssmo > 0,'selected zero ISS model order (non-stationary data or other badness?)');
if ssmo >= ssmomax, fprintf(2,'*** WARNING: selected maximum ISS model order (''ssmomax'' set too low?)\n'); end
fprintf('\nISS model order (SVC) = %d\n',ssmo);

% ISS model estimation

% Estimate ISS model parameters using CCA state-space-subspace algorithm

ptic('\n*** tsdata_to_ss... ');
[A,C,K,V] = tsdata_to_ss(X,pf,ssmo);
ptoc;

% Report information on the estimated ISS, and check for errors

info = ss_info(A,C,K,V);
assert(~info.error,'ISS error(s) found - bailing out');

% Calculate time-domain entropy rate

ERATE = logdet(V);
if ernorm % we subtract the process covariance
	ERATE = ERATE - ss_to_autocov(A,C,K,V,0); % NOTE: this will always be negative - that's fine!
	fprintf('Entropy rate (time domain, normalised) = %g\n',ERATE);
else
	fprintf('Entropy rate (time domain) = %g\n',ERATE);
end

% If not specified, estimate a reasonable frequency resolution for spectral entropy rate

if isempty(fres)
	ptic('\n*** iss2fres... ');
	[fres,frpow2,frierr] = iss2fres(A,C,K,V);
	ptoc;
	fprintf('\nUsing frequency resolution %d = 2^%d (integration error = %.2e)\n',fres,frpow2,frierr);
end

% Calculate cross-power spectral density (CPSD)

S = ss_to_cpsd(A,C,K,V,fres);

return

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
plot_gc({FF,F},{'PWCGC (actual)','PWCGC (estimated)'},[],[maxF maxF]);

%% Granger causality estimation: frequency domain

% Calculate spectral pairwise-conditional causalities resolution from ISS model
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

plot_sgc({ff,f},freqs,'Spectral Granger causalities (blue = actual, red = estimated)');

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
