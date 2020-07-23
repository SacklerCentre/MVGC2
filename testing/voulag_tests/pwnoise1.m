function [y,p] = pwnoise1(n,fs,mtheta,wtheta,gain)

m = n;
%m = 3*n;

x = randn(m,1);
size(x)
fftx = fft(x);

stheta = wtheta*mtheta;
ga = (mtheta^2)/(stheta^2); % Gamma shape parameter
gb = (stheta^2)/mtheta;     % Gamma scale parameter

m1 = length(x)-1;
F = (0:m1)'*(fs/m1);

p = gampdf(F,ga,gb);       % Gamma-distributed frequency peak
p = gain*(std(fftx)/max(p))*p;

y = real(ifft(fftx+p));

%y = y(n+1:2*n);

n
size(y)
