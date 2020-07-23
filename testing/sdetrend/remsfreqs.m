function x = remsfreqs(x,fs,g,h,pad)

if isempty(fs), fs = 2*pi; end % if no sampling frequency supplied assume angular frequency.

assert(isvector(x)  && isnumeric(x),                          'Time series must be a numeric vector or matrix');
assert(isscalar(fs) && isnumeric(fs) && fs > 0,               'Sampling frequency must be a positive scalar')
%assert(isvector(f)  && isnumeric(f)  && all(f > 0 & f < fs/2),'Removal frequencies must be vector of values between zero and the Nyqvist frequency = %g',fs/2);

no = length(x);

xisrv = isrow(x);
if xisrv           % convert input to column vector
	x = x';
end
mu = mean(x);             % save temporal means
x  = bsxfun(@minus,x,mu); % remove means

x = [x;zeros(pad,1)]; % pad
n = length(x);

kg1 = floor(1+n*g(1)/fs);
kg2 = ceil (1+n*g(2)/fs);
kh1 = floor(1+n*(g(1)-h)/fs);
kh2 = ceil (1+n*(g(2)+h)/fs);

sgap = kg2-kg1+1;
sav1 = kg1-kh1;
sav2 = kh2-kg2;

fprintf('gap samples = %3d\n',sgap);
fprintf('av1 samples = %3d\n',sav1);
fprintf('av2 samples = %3d\n\n',sav2);

xft = fft(x);
xfta1 = mean(xft(kh1:(kg1-1)));
xfta2 = mean(xft((kg2+1):kh2));

xft(kg1:kg2) = (xfta1+xfta2)/2;

x = real(ifft(xft));

x = x(1:no);

x = bsxfun(@plus,x,mu); % restore temporal means

if xisrv           % convert output back to row vector
	x = x';
end
