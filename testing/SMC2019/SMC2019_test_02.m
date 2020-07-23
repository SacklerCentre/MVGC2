%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supply: m, LRstats, S

if ~exist('alpha','var'), alpha = 0.05; end
if ~exist('useed','var'), useed = 0;    end

gpterm = 'epsl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SMC2019_testvar5;

if LRstats
	stat = 'LR';
	ptitle = 'Likelihood-ratio statistics\n\n';
else
	stat = 'F';
	ptitle = 'F-statistics\n\n';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rstate = rng_seed(useed);

F  = nan(p,n,n);
Fs  = nan(p,n,n,S);
for s = 1:S
	fprintf('trial %4d of %4d\n',s,S);
	U = varfima_to_tsdata(A,[],[],V,m);
	[~,stats]  = tsdata_to_pwslgc(U,p);
	Fs(:,:,:,s) = significance(1-stats.(stat).cdf(stats.(stat).tstat),alpha,'Bonferroni');
end

rng_restore(rstate);

for i = 1:ncons
	con(i,5) = A(con(i,1),con(i,2),con(i,3));
end
conest  = [];
for i = 1:n
	for j = 1:n
		if j == i, continue; end
		for u = 1:p
			if Fs(u,i,j), conest = [conest;[i,j,u]]; end
		end
	end
end

diag(A(:,:,1))'
con
conest

typeI  = zeros(p,n,n);
typeII = zeros(p,n,n);
for i = 1:n
	for j = 1:n
		if j == i, continue; end
		for u = 1:p
			if A(i,j,u) == 0 % no GC
				typeI (u,i,j) = nnz( Fs(u,i,j,:)); % false positive (Type I)
			else             % yes GC
				typeII(u,i,j) = nnz(~Fs(u,i,j,:)); % false negative
			end
		end
	end
end
