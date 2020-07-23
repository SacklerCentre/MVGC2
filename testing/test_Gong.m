n    = 5;
rho  = 0.9;
w    = 1.0;
g    = 0.5;

%m    = 1000;

%k    = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = specnorm(randn(n),rho);
%A = var_rand(n,p,rho,w);
V = corr_rand(n,g);

X = varfima_to_tsdata(A,[],[],V,m);

Xk = downsample(X,k);

[AA,VV] = tsdata_to_var(X,1,'LWR');

Ak = eye(n);
for i = 1:k
	Ak = A*Ak;
end
