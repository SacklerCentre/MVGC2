% "Harmonic detrend" - remove a specified frequency and harmonics
%
% Inputs
%
% x    Time series (vector or matrix)
% fs   Sample frequency (optional: defaults to angular frequency)
% W    Frequency to remove (in same units as sampling frequency)
% n    Number of harmonics to remove (optional: defaults to ALL)
%
% Output
%
% y    Time series with frequency W and harmonics removed
%
% NOTE 1: Operates column-wise on matrices
%
% NOTE 2: If the number n of harmonics to remove is omitted (or empty) then ALL
% harmonics up to the Nyqvist frequency are removed.

function y = hdetrend1(x,fs,W,n)

if isempty(fs), fs = 2*pi; end                        % if no sampling frequency supplied assume angular frequency.
if nargin < 4 || isempty(n), n = floor(fs/(2*W)); end % if no number of harmonics supplied, remove all up to Nyqvist frequency

assert(ismatrix(x)  && isnumeric(x),                   'Time series must be a numeric vector or matrix');
assert(isscalar(fs) && isnumeric(fs) && fs > 0,        'Bad sampling frequency!')
assert(isscalar(W)  && isnumeric(W),                   'Bad removal frequency!');
assert(isscalar(n)  && isnumeric(n)  && n == floor(n), 'Bad number of harmonics!');
assert(W >0 && W < fs/2,                               'Removal frequency must lie strictly between zero and the Nyqvist frequency = %g',fs/2);
assert(n*W <= fs/2,                                    'Too many harmonics!');

xisrv = isrow(x);

if xisrv % convert input to column vector
	x = x';
end

[T,c] = size(x);   % time series number of time steps, number of columns
mu    = mean(x);   % save temporal means
x     = x-mu;      % remove means
w     = 2*pi*W/fs; % convert to angular frequency

% Build detrend kernel
%
% K(t) = (1/T) sum_{k = 1}^{T-1} cos(k*w*t)
%
% for t = -(T-1) to +(T-1)

wt = w*(-(T-1):(T-1))';
K = ((cos(n*wt) - cos((n+1)*wt))./(1-cos(wt)) - 1)/T;
K(isnan(K)) = 2*n/T; % handle zeros of 1-cos(wt)

% Convolve time series with kernel to obtain detrend compound sinusoid:

s = zeros(size(x));
for i = 1:c
	s(:,i) = conv(K,x(:,i),'valid');
end

% Restore means

y = x-s+mu;

if xisrv % convert output to row vector
	y = y';
end
