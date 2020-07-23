% Sinusoidal detrend - remove specified frequencies and (optionally) harmonics
%
% Inputs
%
% x      Time series (vector or matrix)
% fs     Sample frequency (if empty, defaults to angular frequency on [0,2*pi])
% f      Vector of frequencies to remove in range (0 .. Nyqvist]
% remh   Flag: remove harmonics
%
% Output
%
% x      Time series with frequency f and harmonics removed
%
% Calculates a least-squares sinusoidal fit to time-series data, and removes the sinusoid.
%
% NOTE:  Operates column-wise on matrices

function x = sdetrendh(x,fs,f,remh,g,h,pad)

if isempty(fs), fs = 2*pi; end % if no sampling frequency supplied assume angular frequency.

assert(ismatrix(x)  && isnumeric(x),                             'Time series must be a numeric vector or matrix');
assert(isscalar(fs) && isnumeric(fs) && fs > 0,                  'Sampling frequency must be a positive scalar')
assert(isvector(f)  && isnumeric(f)  && all(f > eps & f <= fs/2),'Removal frequencies must be vector of values between zero and the Nyqvist frequency = %g',fs/2);


adj = nargin > 4;

if remh % add in harmonics of each frequency
	n = floor(fs./(2*f(:)')); % row vector of numbers of harmonics
	nf = length(f);
	w = [];
	for k = 1:nf
		w = [w (1:n(k))*f(k)];
	end
else
	w = f(:)';
end

w = 2*pi*w/fs;     % so w is a (row) vector of angular frequencies to remove

xisrv = isrow(x);
if xisrv           % convert input to column vector
	x = x';
end

[T,c] = size(x);   % time series number of time steps, number of columns
mu    = mean(x);   % save temporal means
x     = bsxfun(@minus,x,mu); % remove means

% Build sinusoidal detrend kernel: K(t) = (2/T) sum_k cos(w(k)*t) on t = -(T-1) ... +(T-1)

K = 2*sum(cos((-(T-1):(T-1))'*w),2)/T;

% Convolve each input column with kernel to obtain a compound least-squares
% sinusoid, and subtract it

s = zeros(size(x));
for i = 1:c
	s(:,i) = conv(K,x(:,i),'valid');
end

if adj

	Tpad = T+pad;
	kg1 = floor(1+Tpad*g(1)/fs);
	kg2 = ceil (1+Tpad*g(2)/fs);
	kh1 = floor(1+Tpad*(g(1)-h)/fs);
	kh2 = ceil (1+Tpad*(g(2)+h)/fs);

	fprintf('gap samples = %3d\n',  kg2-kg1+1);
	fprintf('av1 samples = %3d\n',  kg1-kh1);
	fprintf('av2 samples = %3d\n\n',kh2-kg2);

	xx = [x;zeros(pad,c)]; % pad

	t = (0:T-1)';
	nf = length(f);
	for i = 1:c
		lxfti = localmeanft(xx(:,i),kg1,kg2,kh1,kh2)/T
		for k = 1:nf
			s(:,i) = s(:,i) + 2*real(lxfti)*cos(w(k)*t) + 2*imag(lxfti)*sin(w(k)*t);
		end
	end

end

x = x-s;

x = bsxfun(@plus,x,mu); % restore means

if xisrv           % convert output back to row vector
	x = x';
end

function lxft = localmeanft(x,kg1,kg2,kh1,kh2)

xft = fft(x);
xfta1 = mean(xft(kh1:(kg1-1)));
xfta2 = mean(xft((kg2+1):kh2));
lxft = (xfta1+xfta2)/2;
