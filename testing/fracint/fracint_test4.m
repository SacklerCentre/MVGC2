n = 10;
p = 11;
q = 17;
g = 0.4;
rhoa = 0.95; wa = 0.3;
rhob = 0.85; wb = 0.5;

d = 0.5;
r = 1000;
arform = true;

m = 10000;

A =  var_rand(n,p,rhoa,wa);
B = -var_rand(n,p,rhob,wa);
V =  corr_rand(n,g);

[L,cholp] = chol(V,'lower');
E = L*randn(n,m);

X = varfima_to_tsdata(A,B,[d r arform],E,[],[],0);

gp_qplot((1:m)',X',[],[],[],[1 Inf])
