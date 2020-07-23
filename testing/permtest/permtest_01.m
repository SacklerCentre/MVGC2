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
if ~exist('ares',  'var'), ares  = 1000;  end
if ~exist('seed',  'var'), seed  = 0;     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;
z = nx+ny+1:n;
r  = [x z];

M  = m-p;          % chi^2 scaling factor = effective number of observations
d  = p*nx*ny;      % chi^2 df and F df1
d2 = nx*(M-p*n)-1; % F df2
K  = d2/d;         % F scaling factor

pvals1 = zeros(S,1);
pvals2 = zeros(S,1);
for s = 1:S
	fprintf('sig: sample %4d of %d\n',s,S)
	G = ones(n,n,p); G(x,y,:) = f; % significant
	V = corr_rand(n,q);
	A = var_rand(G,[],rho,w);
	X = var_to_tsdata(A,V,m);
	[~,Vf] = tsdata_to_var(X,     p,'LWR');      % full regression
	[~,Vr] = tsdata_to_var(X(r,:),p,'LWR');      % reduced regression
	FFT  = trace(Vr(x,x))/trace(Vf(x,x)) - 1;    % F test statistic
	FLR  = logdet(Vr(x,x)) - logdet(Vf(x,x));    % likelihood-ratio test statistic
	pvals1(s) = 1-chi2cdf(M*FLR,d);
	pvals2(s) = 1-fcdf(K*FFT,d,d2);
end
pvals1 = sort(pvals1);
pvals2 = sort(pvals2);

fprintf('\n');

pvaln1 = zeros(S,1);
pvaln2 = zeros(S,1);
for s = 1:S
	fprintf('nsg: sample %4d of %d\n',s,S)
	G = ones(n,n,p); G(x,y,:) = 0; % not significant
	V = corr_rand(n,q);
	A = var_rand(G,[],rho,w);
	X = var_to_tsdata(A,V,m);
	[~,Vf] = tsdata_to_var(X,     p,'LWR');      % full regression
	[~,Vr] = tsdata_to_var(X(r,:),p,'LWR');      % reduced regression
	FFT  = trace(Vr(x,x))/trace(Vf(x,x)) - 1;    % F test statistic
	FLR  = logdet(Vr(x,x)) - logdet(Vf(x,x));    % likelihood-ratio test statistic
	pvaln1(s) = 1-chi2cdf(M*FLR,d);
	pvaln2(s) = 1-fcdf(K*FFT,d,d2);
end
pvaln1 = sort(pvaln1);
pvaln2 = sort(pvaln2);

alpha = linspace(0,1,ares);

fpr1 = [0;mean(pvaln1 <= alpha)';1]; % false positive rates
tpr1 = [0;mean(pvals1 <= alpha)';1]; % true  positive rates

fpr2 = [0;mean(pvaln2 <= alpha)';1]; % false positive rates
tpr2 = [0;mean(pvals2 <= alpha)';1]; % true  positive rates

AUC1 = sum((fpr1(2:end)-fpr1(1:end-1)).*(tpr1(2:end)+tpr1(1:end-1)))/2;
AUC2 = sum((fpr2(2:end)-fpr2(1:end-1)).*(tpr2(2:end)+tpr2(1:end-1)))/2;

samps = (1:S)';

gpstem = fullfile(tempdir,'permtest_01');
gp_write([gpstem '_PVS'],[samps sort(pvals1) sort(pvals2)]);
gp_write([gpstem '_PVN'],[samps sort(pvaln1) sort(pvaln2)]);
gp_write([gpstem '_AUC1'],[fpr1 tpr1]);
gp_write([gpstem '_AUC2'],[fpr2 tpr2]);
gp = gp_open(gpstem,'',[Inf 0.7]);
fprintf(gp,'datfilePVS  = "%s_PVS.dat"\n',gpstem);
fprintf(gp,'datfilePVN  = "%s_PVN.dat"\n',gpstem);
fprintf(gp,'datfileAUC1 = "%s_AUC1.dat"\n',gpstem);
fprintf(gp,'datfileAUC2 = "%s_AUC2.dat"\n',gpstem);
fprintf(gp,'\n');
fprintf(gp,'set multiplot title "x = %d, y = %d, z = %d, p = %d, m = %d, f = %.2f" layout 3,1\n',nx,ny,nz,p,m,f);
fprintf(gp,'set key top left\n');
fprintf(gp,'set ylabel "p-values (sorted)" rot\n');
fprintf(gp,'set xlabel "significant"\n');
fprintf(gp,'plot \\\n');
fprintf(gp,'datfilePVS u 1:2 w lines ls 1 t "LR", \\\n');
fprintf(gp,'datfilePVS u 1:3 w lines ls 2 t "FT"    \n');
fprintf(gp,'set xlabel "non-significant"\n');
fprintf(gp,'unset arrow; set arrow from graph 0,graph 0 to graph 1,graph 1 nohead ls 3\n');
fprintf(gp,'plot \\\n');
fprintf(gp,'datfilePVN u 1:2 w lines ls 1 t "LR", \\\n');
fprintf(gp,'datfilePVN u 1:3 w lines ls 2 t "FT"    \n');
fprintf(gp,'set key bottom right\n');
fprintf(gp,'set xlabel "false positive rate"\n');
fprintf(gp,'set ylabel "true positive rate" rot\n');
fprintf(gp,'unset arrow; set arrow from first 0,first 0 to first 1,first 1 nohead ls 3\n');
fprintf(gp,'plot \\\n');
fprintf(gp,'datfileAUC1 u 1:2 w lines ls 1 t "LR test : AUC = %.4f", \\\n',AUC1);
fprintf(gp,'datfileAUC2 u 1:2 w lines ls 2 t "FT test : AUC = %.4f"    \n',AUC2);
fprintf(gp,'unset multiplot\n',gpstem);
gp_close(gp,gpstem);
