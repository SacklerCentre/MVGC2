
T = 5;
m = 9999;

w = 2*pi;

phi = w/10;

%-------------------------------------------------------------

t = linspace(0,T,m+1)';
x = cos(w*t);

y = real(exp(-i*phi)*hilbert(x)); % phase-shift y

gp_qplot(t,[x y])
