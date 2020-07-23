global have_mvfilter;

n = 10;
p = 15;
q = 20;
A = var_specrad(randn(n,n,p),0.8);
%A = squeeze(var_specrad(randn(1,1,p),0.8));
B = var_specrad(randn(n,n,q),0.9);
%B = squeeze(var_specrad(randn(1,1,q),0.9));
m = 100000;
x = randn(n,m);

have_mvfilter = false;
disp('scripted version')
ptic
ym = mvfilter(B,A,x);
ptoc

have_mvfilter = true;
disp('mex version')
ptic
yx = mvfilter(B,A,x);
ptoc

maxabs(ym-yx)
