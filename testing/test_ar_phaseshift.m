a = 0.8;
b = 0.7;
c = 2;

m = 2^16;

maxmo = 60;

regmode = 'LWR';

phi = 0.1;

%-------------------------------------------------------------

AA = [a c; 0 b]
VV = eye(2)
[u,e] = var_to_tsdata(AA,VV,m);

p = size(AA,3)

[A,V] = tsdata_to_var(u,p,regmode)

x = u(1,:);
y = u(2,:);

[AY,VY] = tsdata_to_var(y,p,regmode)

yps = real(exp(-2i*pi*phi)*hilbert(y)); % phase-shift y

ups = [x;yps];

figure(1); clf
[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(ups,maxmo,regmode);

[APS,VPS] = tsdata_to_var(ups,molrt,regmode)
