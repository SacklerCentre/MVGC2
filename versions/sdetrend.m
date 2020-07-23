% Sinusoidal detrend - remove specified frequencies
%
% Inputs
%
% x      Time series (vector or matrix)
% fs     Sample frequency (if empty, defaults to angular frequency on [0,2*pi])
% f      Vector of frequencies to remove in range (0 .. Nyqvist]
% wfac   detrend weighting factor
%
% Output
%
% x      Time series with frequency f removed
%
% Calculates a least-squares sinusoidal fit to time-series data, and removes the sinusoid.
%
% NOTE:  Operates column-wise on matrices

function x = sdetrend(x,fs,f,wfac)

if isempty(fs), fs = 2*pi; end % if no sampling frequency supplied assume angular frequency

if nargin < 4 || isempty(wfac), wfac = 1; end

assert(ismatrix(x)    && isnumeric(x),                                'Time series must be a numeric vector or matrix');
assert(isscalar(fs)   && isnumeric(fs)    && fs > 0,                  'Sampling frequency must be a positive scalar')
assert(isvector(f)    && isnumeric(f)     && all(f > eps & f <= fs/2),'Removal frequencies must be vector of values between zero and the Nyqvist frequency = %g',fs/2);
assert(isvector(wfac) && isnumeric(wfac)  && all(wfac >= 0),          'Weighting factors must be a vector of non-negative values');

w = (2*pi/fs)*f(:)'; % convert w to row vector of angular frequencies
nfreqs = length(w);

xisrv = isrow(x);
if xisrv           % convert input to column vector
	x = x';
end

[T,c] = size(x);   % time series number of time steps, number of columns
mu    = mean(x);   % save temporal means
x     = bsxfun(@minus,x,mu); % remove means

if isscalar(wfac)
	wfac = ones(size(w));
else
	wfac = wfac(:)';
	nwfac = length(wfac);
	if     nwfac < nfreqs
		wfac = [wfac ones(1,nfreqs-nwfac)]; % pad
	elseif nwfac > nfreqs
		wfac = wfac(1:nfreqs);              % truncate
	end
	assert(length(wfac) == nfreqs,'Weighting factors do not match frequencies');
end
%fs*w/(2*pi)
%w/pi
%wfac = wfac*(1-w/pi);

% Build sinusoidal detrend kernel: K(t) = (2/T) sum_k cos(w(k)*t) on t = -(T-1) ... +(T-1)

K = 2*sum(bsxfun(@times,wfac,cos((-(T-1):(T-1))'*w)),2)/T;

% Convolve each input column with kernel to obtain a compound least-squares
% sinusoid, and subtract it (weighted)

for i = 1:c
	x(:,i) = x(:,i) - conv(K,x(:,i),'valid');
end

x = bsxfun(@plus,x,mu); % restore means

if xisrv           % convert output back to row vector
	x = x';
end
