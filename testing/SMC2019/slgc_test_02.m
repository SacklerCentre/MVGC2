nx   = 3;
ny   = 5;
nz    = 2;
p     = 11;
rho   = 0.9;
w     = 1.0;
g     = 0.5;

K     = [2 3 7];
c     = 1;

m     = 1000;

alpha = 0.05;

mhc   = 'NONE';

S     = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;

NK = 1:p; NK(K) = [];

inferns = nan(2,2,2,S);

S10 = round(S/10);
for s = 1:S
	AA = var_rand(n,p,rho,w,false);
	AA(x,y,K)  = c*AA(x,y,K);
	AA(x,y,NK) = 0;
	AA = specnorm(AA,rho);

	VV = corr_rand(n,g);

	U = varfima_to_tsdata(AA,[],[],VV,m);

	stats = slgc(U,x,y,p);

	sigF  = significance(stats.F.pval, alpha,mhc);
	inferns(1,1,1,s) = mean(sigF(K));     % true pos
	inferns(1,2,1,s) = mean(1-sigF(NK));  % true neg
	inferns(2,1,1,s) = mean(sigF(NK));    % false pos
	inferns(2,2,1,s) = mean(1-sigF(K));   % false neg

	sigLR = significance(stats.LR.pval,alpha,mhc);
	inferns(1,1,2,s) = mean(sigLR(K));    % true pos
	inferns(1,2,2,s) = mean(1-sigLR(NK)); % true neg
	inferns(2,1,2,s) = mean(sigLR(NK));   % false pos
	inferns(2,2,2,s) = mean(1-sigLR(K));  % false neg

	if rem(s,S10) == 0, fprintf('.'); end
end
fprintf('\n\n');

infern = mean(inferns,4);

%         +-------+-------+
%         |  pos  |  neg  |
% +-------+-------+-------+
% + true  |  1,1  |  1,2  |
% + false |  2,1  |  2,2  |
% +-------+-------+-------+

fprintf(' F stats +-------+-------+\n');
fprintf('         |  pos  |  neg  |\n');
fprintf(' +-------+-------+-------+\n');
fprintf(' + true  | %5.3f | %5.3f |\n',infern(1,1,1),infern(1,2,1));
fprintf(' + false | %5.3f | %5.3f |\n',infern(2,1,1),infern(2,2,1));
fprintf(' +-------+-------+-------+\n\n');

fprintf('LR stats +-------+-------+\n');
fprintf('         |  pos  |  neg  |\n');
fprintf(' +-------+-------+-------+\n');
fprintf(' + true  | %5.3f | %5.3f |\n',infern(1,1,2),infern(1,2,2));
fprintf(' + false | %5.3f | %5.3f |\n',infern(2,1,2),infern(2,2,2));
fprintf(' +-------+-------+-------+\n\n');
