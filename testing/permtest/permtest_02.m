%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('nx',    'var'), nx    = 3;     end
if ~exist('ny',    'var'), ny    = 5;     end
if ~exist('nz',    'var'), nz    = 0;     end
if ~exist('p',     'var'), p     = 8;     end
if ~exist('rho',   'var'), rho   = 0.9;   end
if ~exist('w',     'var'), w     = 1.0;   end
if ~exist('f',     'var'), f     = 0.1;   end
if ~exist('q',     'var'), q     = 1.0;   end
if ~exist('m',     'var'), m     = 1000;  end
if ~exist('S',     'var'), S     = 1000;  end
if ~exist('N',     'var'), N     = 100;   end
if ~exist('ares',  'var'), ares  = 1000;  end
if ~exist('seed',  'var'), seed  = 0;     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S10 = round(S/10);

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;
z = nx+ny+1:n;
r  = [x z];

M  = m-p;          % chi^2 scaling factor = effective number of observations
d  = p*nx*ny;      % chi^2 df and F df1
d2 = nx*(M-p*n)-1; % F df2
K  = d2/d;         % F scaling factor

o = decorrlags(rho,m);

tic; fprintf('significant model ');
pvalsa = zeros(S,1);
pvalsp = zeros(S,1);
P = zeros(N,1);
for s = 1:S
   	if rem(s,S10) == 0, fprintf('.'); end
	G = ones(n,n,p); G(x,y,:) = f; % significant
	V = corr_rand(n,q);
	A = var_rand(G,[],rho,w);
	X = var_to_tsdata(A,V,m);
	[~,Vf] = tsdata_to_var(X,     p,'LWR');      % full regression
	[~,Vr] = tsdata_to_var(X(r,:),p,'LWR');      % reduced regression
	Xp = X;
	if LR % likelihood-ratio test
		LRfx = logdet(Vf(x,x));
		LRrx = logdet(Vr(x,x));
		F = LRrx - LRfx;    % likelihood-ratio test statistic
		for i = 1:N
			Xp(y,:) = tsdata_rotate(X(y,:),o);
			[~,Vp]  = tsdata_to_var(Xp,p,'LWR'); % permutation full regression (reduced is same)
			LRpx    = logdet(Vp(x,x));
			P(i)    = LRrx - LRpx;    % likelihood-ratio empirical null
		end
		pvalsa(s) = 1-chi2cdf(M*F,d); % analytic
		pvalsp(s) = mean(F <= P);     % permutation test
	else % F test
		FTfx = trace (Vf(x,x));
		FTrx = trace (Vr(x,x));
		F = FTrx/FTfx - 1;  % F test statistic
		for i = 1:N
			Xp(y,:) = tsdata_rotate(X(y,:),o);
			[~,Vp]  = tsdata_to_var(Xp,p,'LWR'); % permutation full regression (reduced is same)
			FTpx    = trace (Vp(x,x));
			P(i)    = FTrx/FTpx - 1;  % F  empirical null
		end
		pvalsa(s) = 1-fcdf(K*F,d,d2); % analytic
		pvalsp(s) = mean(F <= P);     % permutation test
	end
end
pvalsa = sort(pvalsa);
pvalsp = sort(pvalsp);
fprintf(' %s\n',datestr(toc/(60*60*24),'HH:MM:SS'));

tic; fprintf('non-significant model ');
pvalna = zeros(S,1);
pvalnp = zeros(S,1);
P = zeros(N,1);
for s = 1:S
   	if rem(s,S10) == 0, fprintf('.'); end
	G = ones(n,n,p); G(x,y,:) = 0; % not significant
	V = corr_rand(n,q);
	A = var_rand(G,[],rho,w);
	X = var_to_tsdata(A,V,m);
	[~,Vf] = tsdata_to_var(X,     p,'LWR');      % full regression
	[~,Vr] = tsdata_to_var(X(r,:),p,'LWR');      % reduced regression
	if LR % likelihood-ratio test
		LRfx = logdet(Vf(x,x));
		LRrx = logdet(Vr(x,x));
		F = LRrx - LRfx;    % likelihood-ratio test statistic
		for i = 1:N
			X(y,:) = tsdata_rotate(X(y,:),o);
			[~,Vp] = tsdata_to_var(X,p,'LWR'); % permutation full regression (reduced is same)
			LRpx   = logdet(Vp(x,x));
			P(i) = LRrx - LRpx;    % likelihood-ratio empirical null
		end
		pvalna(s) = 1-chi2cdf(M*F,d); % analytic
		pvalnp(s) = mean(F <= P);     % permutation test
	else % F test
		FTfx = trace (Vf(x,x));
		FTrx = trace (Vr(x,x));
		F = FTrx/FTfx - 1;  % F test statistic
		for i = 1:N
			X(y,:) = tsdata_rotate(X(y,:),o);
			[~,Vp] = tsdata_to_var(X,p,'LWR'); % permutation full regression (reduced is same)
			FTpx   = trace (Vp(x,x));
			P(i) = FTrx/FTpx - 1;  % F  empirical null
		end
		pvalna(s) = 1-fcdf(K*F,d,d2); % analytic
		pvalnp(s) = mean(F <= P);     % permutation test
	end
end
pvalna = sort(pvalna);
pvalnp = sort(pvalnp);
fprintf(' %s\n',datestr(toc/(60*60*24),'HH:MM:SS'));

alpha = linspace(0,1,ares);

fpra = [0;mean(pvalna <= alpha)';1]; % false positive rates
tpra = [0;mean(pvalsa <= alpha)';1]; % true  positive rates
fprp = [0;mean(pvalnp <= alpha)';1]; % false positive rates
tprp = [0;mean(pvalsp <= alpha)';1]; % true  positive rates

AUCa = sum((fpra(2:end)-fpra(1:end-1)).*(tpra(2:end)+tpra(1:end-1)))/2;
AUCp = sum((fprp(2:end)-fprp(1:end-1)).*(tprp(2:end)+tprp(1:end-1)))/2;

samps = (1:S)';

if LR,ttype = 'LR-test'; else, ttype = 'F-test'; end;

gpstem = fullfile(tempdir,'permtest_02');
gp_write([gpstem '_PVL'],[samps pvalsa pvalsp pvalna pvalnp]);
gp_write([gpstem '_AUC'],[fpra tpra fprp tprp]);
gp = gp_open(gpstem,'',[Inf 0.7]);
fprintf(gp,'datfilePVL = "%s_PVL.dat"\n',gpstem);
fprintf(gp,'datfileAUC = "%s_AUC.dat"\n',gpstem);
fprintf(gp,'\n');
fprintf(gp,'set multiplot title "%s: x = %d, y = %d, z = %d, p = %d, m = %d, f = %.2f" layout 3,1\n',ttype,nx,ny,nz,p,m,f);
fprintf(gp,'set key top left\n');
fprintf(gp,'set ylabel "p-values (sorted)" rot\n');
fprintf(gp,'set xlabel "significant"\n');
fprintf(gp,'unset arrow; set arrow from graph 0,graph 0 to graph 1,graph 1 nohead ls 3\n');
fprintf(gp,'plot \\\n');
fprintf(gp,'datfilePVL u 1:2 w lines ls 1 t "analytic", \\\n');
fprintf(gp,'datfilePVL u 1:3 w lines ls 2 t "permtest"    \n');
fprintf(gp,'set xlabel "non-significant"\n');
fprintf(gp,'unset arrow; set arrow from graph 0,graph 0 to graph 1,graph 1 nohead ls 3\n');
fprintf(gp,'plot \\\n');
fprintf(gp,'datfilePVL u 1:4 w lines ls 1 t "analytic", \\\n');
fprintf(gp,'datfilePVL u 1:5 w lines ls 2 t "permtest"    \n');
fprintf(gp,'set key bottom right\n');
fprintf(gp,'set xlabel "false positive rate"\n');
fprintf(gp,'set ylabel "true positive rate" rot\n');
fprintf(gp,'unset arrow; set arrow from first 0,first 0 to first 1,first 1 nohead ls 3\n');
fprintf(gp,'plot \\\n');
fprintf(gp,'datfileAUC u 1:2 w lines ls 1 t "analytic test : AUC = %.4f", \\\n',AUCa);
fprintf(gp,'datfileAUC u 3:4 w lines ls 2 t "permtest test : AUC = %.4f"    \n',AUCp);
fprintf(gp,'unset multiplot\n',gpstem);
gp_close(gp,gpstem);
