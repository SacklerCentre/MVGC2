% supply a, b, c

q = 8;
m = 40;

%-------------------------------------------------------------

n = 2;
x = 1;
y = 2;

VARA = zeros(n,n,q);
VARA(:,:,1) = [a 0;0 b];
VARA(:,:,q) = [0 c;0 0];
V    = eye(2);

var_specrad(VARA)

F1 = var_to_mvgc(VARA,V,x,y)
F = var_to_mvgc_ms(VARA,V,x,y,m)

D = (1+b^2+c^2)/2;
v = D+sqrt(D^2-b^2);
h = b/v;

m1 = m-1;
k = (1:m1)';

rk  = (b.^(k+1)-a.^(k+1))/(b-a);
rk1 = [0;rk(1:m1-1)];
rkq = [zeros(q-1,1);rk(1:m-q)];

bkxx = (rk-h*rk1).^2;
Bkxx = (rk-b*rk1).^2 + (c^2)*rkq.^2;

bkxx = [1;bkxx];
Bkxx = [1;Bkxx];

Fm = log(v) + log(cumsum(bkxx)) - log(cumsum(Bkxx));

gp_qplot((1:m)',[F Fm],[],sprintf('set grid\nset xr [1:%d]',m));
