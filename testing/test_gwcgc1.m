group = {[4 3] [9 8] [6 7 5 1 2]};

N = 5;
m = 500;

rng_seed(423)
%-------------------------------------------------------------

g = length(group);

A0 = var9_test;
[n,~,p] = size(A0);
V0 = corr_rand(n,1);

X = var_to_tsdata(A0,V0,m,N);

[A,V] = tsdata_to_var(X,p);

[F1 stats1] = var_to_gwcgc(A,V,group,'both',X);
F1

F2 = nan(g);
stats = cell(g);
for j = 1:g
	for i = 1:g
		if i == j, continue; end
		[F2(i,j),stats{i,j}] = var_to_mvgc(A,V,group{i},group{j},'both',X);
	end
end
F2
F1-F2

statsf = struct('tstat',nan(g),'pval',nan(g));
stats2.single = struct('F',statsf,'LR',statsf);
stats2.dual = struct('F',statsf,'LR',statsf);
clear statsf;
for j = 1:g
	for i = 1:g
		if i == j, continue; end
		stats2.single.F.tstat(i,j)  = stats{i,j}.single.F.tstat;
		stats2.single.LR.tstat(i,j) = stats{i,j}.single.LR.tstat;
		stats2.single.F.pval(i,j)   = stats{i,j}.single.F.pval;
		stats2.single.LR.pval(i,j)  = stats{i,j}.single.LR.pval;
		stats2.dual.F.tstat(i,j)    = stats{i,j}.dual.F.tstat;
		stats2.dual.LR.tstat(i,j)   = stats{i,j}.dual.LR.tstat;
		stats2.dual.F.pval(i,j)     = stats{i,j}.dual.F.pval;
		stats2.dual.LR.pval(i,j)    = stats{i,j}.dual.LR.pval;
	end
end
clear stats

fprintf('\n');
fprintf('|F1-F2| = %g\n',maxabs(F1-F2));
