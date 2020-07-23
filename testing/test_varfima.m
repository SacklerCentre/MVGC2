global have_mvfilter;
have_mvfilter = true;

n = 10;
p = 15;
q = 20;
A = var_specrad(randn(n,n,p),0.8);
%A = squeeze(var_specrad(randn(1,1,p),0.8));
B = var_specrad(randn(n,n,q),0.9);
%B = squeeze(var_specrad(randn(1,1,q),0.9));
m = 100000;
x = randn(n,m);

fiparms = [0.4,1000,true];

ptic
y = varfima_to_tsdata(A,B,fiparms,x);
ptoc

ptic
y1 = varfima_to_tsdata1(A,B,fiparms,x);
ptoc

maxabs(y1-y)
