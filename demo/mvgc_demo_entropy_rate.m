%%%%%%%%%%%%%%%%%%%%%%%%% State-space entropy rate demo %%%%%%%%%%%%%%%%%%%%%%%%

% This script demonstrates how to calculate the entropy rate of a time series
% via state-space modelling, as described in:
%
%   Spectrally and temporally resolved estimation of neural signal diversity (2023),
%   Pedro A.M. Mediano, Fernando E. Rosas, Andrea I. Luppi, Valdas Noreika,
%   Anil K. Seth, Robin L. Carhart-Harris, Lionel Barnett and Daniel Bor,
%   bioRxiv 2023.03.30.534922; doi: https://doi.org/10.1101/2023.03.30.534922
%   (In review, eLife, Nov. 2023.)
%
% Entropy rates are calculated under Gaussian assumptions in the time and frequency
% (spectral) domains, and "band-limited" (i.e., for specified frequency bands), via
% state-space modelling. The units are nats; since entropy rates depend on the
% sampling frequency, it arguably makes more sense to convert results into
% nats-per-second; to achieve this, the entropy rates calculated here should be
% multiplied by the sampling rate (fs).
%
% For a stationary process, entropy rate measures the entropy (unpredictability) of
% the process conditional on its (infinite) past. Note that for continuous-state data
% these are *differential* entropies, and may therefore be negative, which somewhat
% compromises the interpretation as "unpredictability". Differential entropy is also
% not scale-independent; multiplying a time series by a constant factor changes the
% entropy rate. To address this situation, we offer the option to calculate instead
% a "normalised negentropy rate"; in the time domain, this is the mutual information
% between the past and (1-step) future of the process - the *predictability* rather
% than unpredictability. It has the advantage (at least in the time domain) of being
% non-negative and scale-independent, and thereby more clearly interpretable.
%
% Note that entropy rate for a multivariate time series is not calculated per-variable;
% it is an "overall" entropy rate. If you want per-variable (channel) entropy rates,
% you must run this routine seperately for each individual channel.
%
% You may use this script as a template to calculate entropy rates for your data; by
% default, synthetic test data is generated (note that state-space model estimation
% may potentially fail for the generated data - or indeed for your own data). To
% calculate entropy rates for your own multivariate, possibly epoched time-series data
% see comments in the code below.
%
%%%%%%%%%%%%%%%%%%%%%%%%% Parameters (override on command line) %%%%%%%%%%%%%%%%

% Test data generation parameters [NOTE: omit if supplying your own time series data]

defvar('nchans',    10    ); % number of channels
defvar('nepochs',   20    ); % number of epochs (may be 1)
defvar('nobs',      1000  ); % number of observations per epoch
defvar('ssmoact',   25    ); % ISS model order
defvar('rho',       0.95  ); % AR spectral radius
defvar('rmi',       1     ); % residuals log-generalised correlation (multi-information)
defvar('seed',      0     ); % random seed (0 for unseeded)

% Other parameters

defvar('fs',        200   ); % sample rate (Hz)
defvar('nne',       false ); % calculate "normalised negentropy rate" rather than entropy rate.
defvar('varmomax',  32    ); % maximum model order for VAR model order selection (required for ISS estimation)
defvar('fres',      []    ); % frequency resolution  for spectral entropy rate (empty for automatic calculation)

%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN GENERATE TEST DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random ISS test data: FOR YOUR OWN DATA, OMIT AND SKIP TO SECTION BELOW

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

%%%%%%%%%%%%%%%%%%%%%%%%%% END GENERATE TEST DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN LOAD MY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% X = load_my_data(...); % or whatever... dimensions: [channels x observations x epochs]
% [nchans,nobs,nepochs] = size(X);
%
% If your data is not epoched, dimensions of X should be just [channels x observations].
%
%%%%%%%%%%%%%%%%%%%%%%%%%% END LOAD MY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VAR model order estimation (required for ISS estimation)

ptic('\n*** tsdata_to_varmo... ');
varmo = tsdata_to_varmo(X,varmomax,'LWR',[],false,1); % returns optimal model order according to the AIC
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
if ssmo >= ssmomax, fprintf(2,'*** WARNING: selected maximum ISS model order (''pf'' set too low?)\n'); end
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

erate_td = logdet(V)/2;
if nne % normalised negentropy rate: subtract entropy rate from process entropy
	eproc    = logdet(ss_to_autocov(A,C,K,V,0))/2; % process entropy
	erate_td = eproc - erate_td;                   % NOTE: this will be non-negative
	fprintf('Time-domain normalised negentropy rate = %g\n',erate_td);
else
	fprintf('Time-domain entropy rate = %g\n',erate_td);
end

% If not specified, estimate a reasonable frequency resolution for spectral entropy rate

if isempty(fres)
	ptic('\n*** iss2fres... ');
	[fres,frpow2,frierr] = iss2fres(A,C,K,V);
	ptoc;
	fprintf('\nUsing frequency resolution %d = 2^%d (integration error = %.2e)\n',fres,frpow2,frierr);
end

% Calculate cross-power spectral density (CPSD) to Nyqvist frequency at specified resolution

S = ss_to_cpsd(A,C,K,V,fres);

% Calculate frequency-domain (spectral) entropy rate, and plot it

ptic('\n*** calculating spectral entropy rate... ');
erate_fd = zeros(fres+1,1);
for k = 1:fres+1
	erate_fd(k) = logdet(S(:,:,k))/2;
end
ptoc;

if nne % normalised negentropy rate: subtract entropy rate from process entropy
	erate_fd = eproc-erate_fd; % NOTE: this will not necessarily be non-negative!
end

% Sanity check: spectral entropy rate should average (approximately) to time-domain value

abserrs = abs(erate_td-bandlimit(erate_fd,[],fs));
fprintf('\nSpectral integral check: abs. error = %.2e\n',abserrs);
if abserrs > sqrt(eps)
	fprintf(2,'*** WARNING: error is a bit large\n');
end

figure(3); clf;
f = sfreqs(fres,fs); % get frequency vector
plot(f,erate_fd);
if nne
	title('Spectral normalised negentropy rate');
else
	title('Spectral entropy rate');
end
xlabel('Frequency (Hz)');
ylabel('entropy rate (nats)');
xlim([0,fs/2]); % zero to Nyqvist frequency
yline(erate_td,'r'); % red line at time-domain value
grid on

% Band-limited entropy rates (standard frequency bands)

fbands = [0,4;      ... % delta
          4,8;      ... % theta
          8,15;     ... % alpha
          15,30;    ... % beta
          30,50;    ... % low gamma
          50,fs/2];     % high gamma

ptic('\n*** calculating band-limited entropy rates... ');
nfbands = size(fbands,1);
erate_bl = zeros(nfbands,1);
for i = 1:nfbands
	erate_bl(i) = bandlimit(erate_fd,[],fs,fbands(i,:));
end
ptoc;

if nne
	fprintf('\nBand-limited normalised negentropy rates (nats)\n');
	fprintf('------------------------------------------------\n');
else
	fprintf('\nBand-limited entropy rates (nats)\n');
	fprintf('-------------------------------------\n');
end
fprintf('delta      : % 6.4f\n',erate_bl(1));
fprintf('theta      : % 6.4f\n',erate_bl(2));
fprintf('alpha      : % 6.4f\n',erate_bl(3));
fprintf('beta       : % 6.4f\n',erate_bl(4));
fprintf('low gamma  : % 6.4f\n',erate_bl(5));
fprintf('high gamma : % 6.4f\n',erate_bl(6));
if nne
	fprintf('------------------------------------------------\n');
else
	fprintf('-------------------------------------\n');
end

% Sanity check: band-limited entropy rates should sum (approximately) to time-domain value

abserrb = abs(erate_td-sum(erate_bl));
fprintf('\nBand-limited sum check: abs. error = %.2e\n\n',abserrb);
if abserrb > sqrt(eps)
	fprintf(2,'*** WARNING: error is a bit large\n');
end
