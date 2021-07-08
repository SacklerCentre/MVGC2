% SINUSOIDAL DETREND
%
% INVOCATION
%
% [y,s] = sinudet(x,fs,f,sw)
%
% INPUTS
%
% x         matrix of time-series values (variables x observations)
% fs        sampling frequency (default: angular frequency in range 0 ... 2*pi)
% f         sinusoidal frequency to remove (in range (0 .. Nyqvist])
% sw        weights on sinusoids (default: 1 - remove entire sinusoid)
%
% OUTPUTS
%
% y         the sinusoidally-detrended signal
% s         the sinusoidal least-squares fit signal
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [y,s] = sinudet(x,fs,f,sw)

if isempty(fs), fs = 2*pi; end % If no sampling frequency supplied assume angular frequency.

if nargin < 4 || isempty(sw), sw = 1; end % default: remove entire sinusoid

assert(ismatrix(x),'Time series must be a matrix');
[n,m] = size(x);

assert(isvector(f),'Sinusoidal frequencies must be a vector');
assert(all(f > 0 & f <= fs/2),'Sinusoidal frequencies must lie in range (0 .. Nyqvist]');
nf = size(f,1);

if isscalar(sw)
	sw = sw*ones(n,nf);
elseif isrow(sw)
	assert(length(sw) == nf,'Sinusoidal weights must be conformant with time series');
	sw = repmat(sw,n,1);
elseif iscol(sw)
	assert(length(sw) == n,'Sinusoidal weights must be conformant with time series');
	sw = repmat(sw,1,nf);
elseif ismatrix(sw)
	assert(isequal(size(sw),[n nf]),'Sinusoidal weights must be conformant with time series');
else
	error('Sinusoidal weights must be a scalar, vector, or matrix');
end

mu = mean(x,2); % save mean
x  = x-mu;      % temporal de-mean
t  = 0:m-1;     % time sequence
w  = 2*pi*f/fs; % convert to angular frequency

% Weighted sum of sinusoids

s = zeros(size(x));
for k = 1:nf
	s = s + sw(:,k).*sinusig(w(k));
end

y = x-s+mu; % the sinusoidally-detrended x (reinstate mean)

% Nested function to calculate sinusoidal signal

    function ss = sinusig(ww)

	W    = (sin(ww*m)/sin(ww))/m;
	u    = W*cos(ww*(m-1));
	v    = W*sin(ww*(m-1));
	comt = cos(ww*t);
	somt = sin(ww*t);
	xc   = mean(x.*comt,2);
	xs   = mean(x.*somt,2);
	D    = 2/(1-W*W);
	p    = D*((1-u)*xc-v*xs); % cos coefficient
	q    = D*((1+u)*xs-v*xc); % sin coefficient
	ss   = p.*comt+q.*somt;   % the sinusoidal signal

	end % function sinusig

end % function sinudet
