% Sinusoidal detrend - remove specified frequencies and (optionally) harmonics
%
% Inputs
%
% x      Time series (vector or matrix)
% fs     Sample frequency (if empty, defaults to angular frequency on [0,pi])
% f      Vector of frequencies to remove (in same units as sampling frequency)
% remh   Flag: remove harmonics
% wfac   detrend weighting factor
%
% Output
%
% x      Time series with frequency f and harmonics removed
%
% Calculates a least-squares sinusoidal fit to time-series data, and removes the sinusoid.
%
% NOTE:  Operates column-wise on matrices

function x = sdetrend_with_harmonics(x,fs,f,remh,wfac)

if isempty(fs), fs = 2*pi; end % if no sampling frequency supplied assume angular frequency.

assert(ismatrix(x)  && isnumeric(x),                          'Time series must be a numeric vector or matrix');
assert(isscalar(fs) && isnumeric(fs) && fs > 0,               'Sampling frequency must be a positive scalar')
assert(isvector(f)  && isnumeric(f)  && all(f > 0 & f < fs/2),'Removal frequencies must be vector of values between zero and the Nyqvist frequency = %g',fs/2);

if remh % add in harmonics of each frequency (up to Nyqvist)
	n = floor(fs./(2*f(:)')); % row vector of numbers of harmonics
	nf = length(f);
	w = [];
	for k = 1:nf
		w = [w (1:n(k))*f(k)];
	end
else
	w = f(:)';
end

if nargin < 5 || isempty(wfac), wfac = 1; end

w = 2*pi*w/fs;     % so w is a (row) vector of angular frequencies to remove

xisrv = isrow(x);
if xisrv           % convert input to column vector
	x = x';
end

[T,c] = size(x);   % time series number of time steps, number of columns
mu    = mean(x);   % save temporal means
x     = bsxfun(@minus,x,mu); % remove means

%fs*w/(2*pi)
%w/pi
wfac = wfac*(1-w/pi);


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
