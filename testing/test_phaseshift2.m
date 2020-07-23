
T = 5;
m = 1000;

a = 0.99;

phi = (2*pi)/10;

%-------------------------------------------------------------

t = linspace(0,T,m)';

x = var_to_tsdata(a,1,m)';

y = real(exp(-i*phi)*hilbert(x)); % phase-shift y

gp_qplot(t,[x y])
