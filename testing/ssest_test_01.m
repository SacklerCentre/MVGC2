%% MVGC state-space demo

%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test data generation

ntrials   = 10;     % number of trials
nobs      = 10000;  % number of observations per trial
fs        = 200;    % sample rate (Hz)

% Actual VAR model generation parameters

nvars     = 5;      % number of variables
ssmoact   = 11;      % SS model order
rhoa      = 0.95;   % AR spectral radius
wvar      = 0.5;    % var coefficients decay weighting factor
rmi       = 0.5;    % residuals log-generalised correlation (multi-information)
                    % g = -log|R|. g = 0 yields zero correlation,g = [] is uniform random
                    % on space of correlation matrices

% VAR model order estimation

varmosel  = 'LRT';  % VAR model order selection ('ACT', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
varmomax  = nvars*ssmoact; % maximum model order for VAR model order selection

% SS model order estimation

ssmosel   = 'SVC';  % SS model order selection ('ACT', 'SVC', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)

% MVGC (frequency domain)

fres      = [];     % spectral MVGC frequency resolution (empty for automatic calculation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('seed',   'var'), seed     = 0;    end % random seed (0 for unseeded)
if ~exist('svconly','var'), svconly  = true; end % only compute SVC for SS model order selection (faster)
if ~exist('plotm',  'var'), plotm    = 0;    end % plot mode (figure number offset, or Gnuplot terminal string)

%% Generate random SS test data

% Seed random number generator.

rng_seed(seed);

% Generate random SS parameters in innovations form

[AA,CC,KK] = iss_rand(nvars,ssmoact,rhoa);

% Generate random residuals covariance (in fact correlation) matrix.

VV = corr_rand(nvars,rmi);

% Report information on the generated VAR model and check for errors.

infoo = ss_info(AA,CC,KK,VV);
assert(~infoo.error,'SS error(s) found - bailing out');

% Generate multi-trial SS time series data with normally distributed residuals
% for generated SS parameters and residuals covariance matrix.

ptic('*** ss_to_tsdata... ');
X = ss_to_tsdata(AA,CC,KK,VV,nobs,ntrials);
ptoc;

%% VAR model order estimation

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.

ptic('\n*** tsdata_to_varmo... ');
if isnumeric(plotm), plotm = plotm+1; end
[varmoaic,varmobic,varmohqc,varmolrt] = tsdata_to_varmo(X,varmomax,'LWR',[],[],plotm);
ptoc;

% Select and report VAR model order.

varmo = moselect(sprintf('VAR model order selection (max = %d)',varmomax),varmosel,'AIC',varmoaic,'BIC',varmobic,'HQC',varmohqc,'LRT',varmolrt);
assert(varmo > 0,'selected zero model order! GCs will all be zero!');
if varmo >= varmomax, fprintf(2,'*** WARNING: selected VAR maximum model order (may have been set too low)\n'); end

%% SS model order estimation

pf = 2*varmo; % Bauer recommends 2 x VAR AIC model order

if svconly  % SVC only: computationally much faster

	ptic('\n*** tsdata_to_sssvc... ');
	if isnumeric(plotm), plotm = plotm+1; end
	[ssmosvc,ssmomax] = tsdata_to_sssvc(X,pf,[],plotm);
	ptoc;

	% Select and report SS model order.

	ssmo = moselect(sprintf('SS model order selection (max = %d)',ssmomax),ssmosel,'ACT',ssmoact,'SVC',ssmosvc);


else        % SVC + likelihood-based selection criteria + SVC: computational intensive

	ptic('\n*** tsdata_to_ssmo... ');
	if isnumeric(plotm), plotm = plotm+1; end
	[ssmoaic,ssmobic,ssmohqc,ssmosvc,ssmolrt,ssmomax] = tsdata_to_ssmo(X,pf,[],[],plotm);
	ptoc;

	% Select and report SS model order.

	ssmo = moselect(sprintf('SS model order selection (max = %d)',ssmomax),ssmosel,'ACT',ssmoact,'AIC',ssmoaic,'BIC',ssmobic,'HQC',ssmohqc,'SVC',ssmosvc,'LRT',ssmolrt);

end

assert(ssmo > 0,'selected zero model order! GCs will all be zero!');
if ssmo >= ssmomax, fprintf(2,'*** WARNING: selected SS maximum model order (may have been set too low)\n'); end

%% SS model estimation

% Estimate SS model order and model paramaters

[A,C,K,V] = tsdata_to_ss(X,pf,ssmo);

% Report information on the estimated SS, and check for errors.

info = ss_info(A,C,K,V);
assert(~info.error,'SS error(s) found - bailing out');

if isempty(fres)
    fres = 2^nextpow2(max(info.acdec,infoo.acdec)); % alternatively, fres = 2^nextpow2(nobs);
	fprintf('\nUsing frequency resolution %d\n',fres);
end
if fres > 10000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

S(:,:,:,1) = ss_to_cpsd(AA,CC,KK,VV,fres);
S(:,:,:,2) = ss_to_cpsd(A ,C ,K ,V ,fres);
if isnumeric(plotm), plotm = plotm+1; figure(plotm); clf; end
plot_cpsd(S,{'act','est'},fs,[],false);

h = fres+1;
LDS1 = zeros(h,1);
LDS2 = zeros(h,1);
for k = 1:h; LDS1(k) = logdet(S(:,:,k,1)); end
for k = 1:h; LDS2(k) = logdet(S(:,:,k,2)); end
lam = sfreqs(fres,fs)';
if isnumeric(plotm), plotm = plotm+1; figure(plotm); clf; end
plot(lam,[LDS1 LDS2]);
legend({'act','est'});
title(['\Delta(log|V|) = ' sprintf('%g',logdet(VV)-logdet(V))]);
