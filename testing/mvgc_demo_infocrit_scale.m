
% Test data generation

ntrials   = 4;      % number of trials
nobs      = 1000;   % number of observations per trial
fs        = 200;    % sample rate (Hz)
seed      = 0;      % random seed (0 for unseeded)

% Actual VAR model generation parameters

tnet      = tnet5;  % connectivity network
moact     = 6;      % model order
rho       = 0.95;   % spectral radius
wvar      = 0.5;    % var coefficients decay weighting factor
rmi       = 0.5;    % residuals log-generalised correlation (multi-information)
                    % g = -log|R|. g = 0 yields zero correlation,g = [] is uniform random
                    % on space of correlation matrices

% VAR model order estimation

momax     = 2*moact; % maximum model order for model order selection
moalpha   = 0.01;   % significance level for model order selection tests
motstat   = '';     % statistical test for LR model order selection:  'F' for F-test or 'chi2' for chi2 test (default)

% Seed random number generator.

rng_seed(seed);

% Generate random VAR coefficients for test network.

AA = var_rand(tnet,moact,rho,wvar);
nvars = size(AA,1); % number of variables

% Generate random residuals covariance (in fact correlation) matrix.

VV = corr_rand(nvars,rmi);

% Report information on the generated VAR model and check for errors.

infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

% Generate multi-trial VAR time series data with normally distributed residuals
% for generated VAR coefficients and residuals covariance matrix.

ptic('*** var_to_tsdata... ');
X = var_to_tsdata(AA,VV,nobs,ntrials);
ptoc;
X1 = X;
X2 = a*X+mu;
clear X

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.

figure(1); clf;
ptic('\n*** tsdata_to_varmo... ');
[moaic1,mobic1,mohqc1,molrt1] = tsdata_to_varmo(X1,momax,'OLS',moalpha,motstat);
ptoc;

figure(2); clf;
ptic('\n*** tsdata_to_varmo... ');
[moaic2,mobic2,mohqc2,molrt2] = tsdata_to_varmo(X2,momax,'OLS',moalpha,motstat);
ptoc;
