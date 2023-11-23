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

% If not specified, estimate a reasonable frequency resolution for spectral entropy rate

if isempty(fres)
	ptic('*** iss2fres... ');
	[fres,frpow2,frierr] = iss2fres(A,C,K,V);
	ptoc;
	fprintf('\nUsing frequency resolution %d = 2^%d (integration error = %.2e)\n',fres,frpow2,frierr);
end

% Calculate time-domain entropy rate

ERATE = logdet(V);
if ernorm % we subtract the process covariance
	EPROC = logdet(ss_to_autocov(A,C,K,V,0)); % process entropy
	ERATE = ERATE - EPROC;                    % NOTE: this will always be negative - that's fine!
	fprintf('\nTime-domain entropy rate (nats/sec, normalised) = %g\n',ERATE);
else
	fprintf('\nTime-domain entropy rate (nats/sec) = %g\n',ERATE);
end

% Calculate cross-power spectral density (CPSD)

S = ss_to_cpsd(A,C,K,V,fres);

% Calculate spectral entropy rate, and plot it

ptic('\n*** calculating spectral entropy rate... ');
erate = zeros(fres+1,1);
for k = 1:fres+1
	erate(k) = logdet(S(:,:,k));
end
ptoc;

if ernorm % we subtract the process covariance
	erate = erate-EPROC; % NOTE: this will NOT necessarily be negative - that's fine!
end

% Sanity check: spectral entropy rate should average (approximately) to time-domain value

fprintf('\nSpectral integral check: rel. error = %.2e\n',abs(1-bandlimit(erate,[],fs)/ERATE));

figure(3); clf;
f = sfreqs(fres,fs); % get frequency vector
plot(f, erate);
if ernorm
	title('Spectral entropy rate (normalised)');
else
	title('Spectral entropy rate');
end
xlabel('Frequency (Hz)');
ylabel('entropy rate (nats/sec)');
xlim([0,fs/2]); % zero to Nyqvist frequency
yline(ERATE,'r'); % red line at time-domain value
grid on

% band-limited entropy rates (standard frequency bands)

fbands = [0,4;      ... % delta
          4,8;      ... % theta
          8,15;     ... % alpha
          15,30;    ... % beta
          30,50;    ... % low gamma
          50,fs/2];     % high gamma

ptic('\n*** calculating band-limited entropy rates... ');
nfbands = size(fbands,1);
BLERATE = zeros(nfbands,1);
for i = 1:nfbands
	BLERATE(i) = bandlimit(erate,[],fs,fbands(i,:));
end
ptoc;

if ernorm
	fprintf('\nBand-limited entropy rates (nats/sec,normalised)\n');
	fprintf('------------------------------------------------\n');
else
	fprintf('\nBand-limited entropy rates (nats/sec)\n');
	fprintf('-------------------------------------\n');
end
fprintf('delta      : % 6.4f\n',BLERATE(1));
fprintf('theta      : % 6.4f\n',BLERATE(2));
fprintf('alpha      : % 6.4f\n',BLERATE(3));
fprintf('beta       : % 6.4f\n',BLERATE(4));
fprintf('low gamma  : % 6.4f\n',BLERATE(5));
fprintf('high gamma : % 6.4f\n',BLERATE(6));
if ernorm
	fprintf('------------------------------------------------\n');
else
	fprintf('-------------------------------------\n');
end

% Sanity check: band-limited entropy rates should sum (approximately) to time-domain value

fprintf('\nBand-limited sum check: rel. error = %.2e\n\n',abs(1-sum(BLERATE)/ERATE));
