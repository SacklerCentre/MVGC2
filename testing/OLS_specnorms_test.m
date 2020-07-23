
n = 5;
p = 7;
w = 0.5;
r = 0.5;

m = 600;

rho = 0.99;

S = 100;

%-------------------------------------------------------------------------------

pmax = 3*p;
ps = (1:pmax)';

S10 = round(S/10);
for s = 1:S
   	if rem(s,S10) == 0, fprintf('.'); end % progress indicator

	A = var_rand(n,p,rho,w);
	V = corr_rand(n,r);

	X = var_to_tsdata(A,V,m);

	rhos(:,s) = OLS_specnorms(X,pmax);
end

rhomean = mean(rhos,2);
rhosdev = std(rhos,[],2);

linep = sprintf('set arrow from first %d,graph 0 to first %d,graph 1 nohead',p,p);
line1 = sprintf('set arrow from graph 0,first 1 to graph 1,first 1 nohead');

gp_qplot(ps,[rhomean rhomean+rhosdev rhomean-rhosdev],[],sprintf('unset key\nset grid\n%s\n%s',linep,line1));
