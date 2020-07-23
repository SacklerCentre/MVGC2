% Welch algorithm, with sensible defaults

function [S,f,nwobs,noobs,nwins] = tsdata_to_cpsd(X,fs,window,overlap,fres,autospec)

% Defaults

def_nwobs = 128;
def_nwins = 8;
def_folap = 0.5;
def_fres  = 1024;

[nvars,nobs,ntrials] = size(X);

if nargin < 2 || isempty(fs), fs = 2*pi; end % default to angular frequency

if nargin < 3 || isempty(window)  % default:  number of observations per window
	have_nwins = false;
	nwobs = def_nwobs;
elseif window > eps               % positive: number of observations per window
	have_nwins = false;
	nwobs = window;
elseif window < -eps              % negative: number of windows (with overlap)
	have_nwins = true;
	nwins = -window;
else                              % zero:     number of windows (with overlap)
	have_nwins = true;
	nwins = def_nwins;
end

if nargin < 4 || isempty(overlap) % default:  fraction of window to overlap
	have_folap = true;
	folap = def_folap;
elseif  overlap > eps             % positive: fraction of window to overlap
	have_folap = true;
	folap = overlap;
elseif  overlap < -eps            % negative: number of observations per overlap
	have_folap = false;
	noobs = -overlap;
else
	error('bad window/overlap!');
end

% Make sure we have nwobs and noobs

if have_nwins
	if have_folap
		noobs = round(folap*round(nobs/(nwins-folap*(nwins-1))));
	end
	nwobs = floor((nobs+(nwins-1)*noobs)/nwins);
else
	if have_folap
		noobs = round(folap*nwobs);
	end
	nwins = floor((nobs-noobs)/(nwobs-noobs));
end
assert(noobs > 0 && noobs < nwobs,'bad window/overlap!');

% Make sure we have nfft

if nargin < 5 || isempty(fres) % pwelch default nfft not necessarily good for spectral factorisation...
	nfft = 2*def_fres;
else
	if fres > eps
		nfft = 2*fres;
	elseif fres < -eps
		error('bad frequency resolution!');
	else                        % (zero) ...but if you really want it, set fres = 0
		nfft = max(256,2^(nextpow2(nwobs)));
	end
end

if nargin < 6 || isempty(autospec), autospec = false; end

fres = nfft/2;
h = fres+1;
if autospec
	S = zeros(h,nvars);
else
	S = zeros(h,nvars,nvars);
end
f = linspace(0,fs/2,h)';

X = permute(X,[2 1 3]);

if autospec
	for r = 1:ntrials
		for i = 1:nvars
			S(:,i) = S(:,i) + pwelch(X(:,i,r),nwobs,noobs,nfft,fs);
		end
	end
	S = (fs/2)*S/ntrials;
	S(1  ,:) = 2*S(1,  :);
	S(end,:) = 2*S(end,:);
else
	for r = 1:ntrials
		for i = 1:nvars
			S(:,:,i) = S(:,:,i) + cpsd(X(:,i,r),X(:,:,r),nwobs,noobs,nfft,fs);
		end
	end
	S = (fs/2)*permute(S/ntrials,[3 2 1]);
	S(:,:,1  ) = 2*S(:,:,1  );
	S(:,:,end) = 2*S(:,:,end);
end
