p   = 3;
rho = 0.9;
w   = 0.5;
g   = 1;

m = 10000;
N = 10;

regmode = 'OLS';

%rng_seed(21312);

%-------------------------------------------------------------

group = {[3 2],[6 9 5 11 10],[4 1 8 7]};
n = max(cell2mat(group));
ng = length(group);

AA = var_rand(n,p,rho,w);
VV = corr_rand(n,g);

X = var_to_tsdata(AA,VV,m,N);

[A,V] = tsdata_to_var(X,p,regmode);

[F,stats] = var_to_gwiogc(A,V,group,'in','both',X,regmode);

FF = nan(ng,1);
for a = 1:ng
	x = group{a};
	y = 1:n; y(x) = [];
%	y = group{a};
%	x = 1:n; x(y) = [];
	[FF(a) statss(a)] = var_to_mvgc(A,V,x,y,'both',X,regmode);
end

[FF F]
