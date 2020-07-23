%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx        = 2;
ny        = 3;
nz        = 1;

r         = 21;     % SS model order
rhoa      = 0.9;    % AR spectral radius
rmi       = 0.3;    % residuals log-generalised correlation
nrmi      = 0.5;    % noise log-generalised correlation (multi-information)

lmin = 0;
lmax = 5;
lres = 100;

if ~exist('seed', 'var'), seed = 0; end % random seed (0 for unseeded)

gpterm = 'epsl';
gpname = fullfile(tempdir,'snrtest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed);

n = nx+ny+nz;
prm = randperm(n);
x = prm(1:nx);
y = prm(nx+1:nx+ny);
z = prm(nx+ny+1:n);

[A,C,K] = iss_rand(n,r,rhoa);
V = corr_rand(n,rmi);
ss_info(A,C,K,V);

F = ss_to_mvgc(A,C,K,V,x,y);

L = chol(V,'lower');
KL = K*L;
V = L*L';
KV = KL*L';
Q = KL*KL';

O = dlyap(A,Q);
W = C*O*C' + V; % cov matrix of reference signal

N = corr_rand(n,nrmi); % cov matrix of noise signal
N = mean(diag(W)./diag(N))*N;

I = cov_to_mvmi(W,x,y); % (x,y|z) partial correlation for ref signal

J = cov_to_mvmi(N,x,y); % (x,y|z) partial correlation for noise signal

lams = linspace(lmin,lmax,lres)';
lams2 = lams.*lams;
FF = zeros(lres,1);
II = zeros(lres,1);
for i = 1:lres
	fprintf('.');
	lam = lams(i);
	lam2 = lams2(i);
	%fprintf('iter %3d of %3d : lambda = %7.4f\n',i,lres,lam);

	CC = lam*C;
	SS = lam*KV;
	RR = lam2*V+N;

	[KK,VV] = ss2iss(A,CC,Q,RR,SS);

	FF(i) = ss_to_mvgc(A,CC,KK,VV,x,y);

	WW = CC*O*CC' + VV; % cov matrix of composite signal

	II(i) = cov_to_mvmi(WW,x,y);
end
fprintf('\n\n');

%SNR = lams2;
SNR = 20 * log10(lams); % dB

hline1 = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead lc "red"\n',F,F);
htext1 = sprintf('set label "pure signal GC" at graph 0.1,first %g left front boxed\n',F);
hline2 = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead lc "red"\n',I,I);
htext2 = sprintf('set label "pure signal MI" at graph 0.1,first %g left front boxed\n',I);
hline3 = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead lc "green"\n',J,J);
htext3 = sprintf('set label "pure noise MI" at graph 0.1,first %g left front boxed\n',J);
cmds = 'unset key\nset xlabel "high noise <--- SNR (dB) ---> low noise"\nunset ylabel\nset style textbox opaque noborder\n';

yr1 = sprintf('set yr [0:%g]\n',1.05*F);
yr2 = sprintf('set yr [%g:%g]\n',0.9*min(min(II),min(I,J)),1.05*max(max(II),max(I,J)));

gp_qplot(SNR,FF,[],['set title "Granger causality"\n'  cmds yr1 hline1 htext1],gpterm,0.8,[],[gpname '_GC']);
gp_qplot(SNR,II,[],['set title "Mutual information"\n' cmds yr2 hline2 htext2 hline3 htext3],gpterm,0.8,[],[gpname '_MI']);
%gp_qplot(lams2,SNR);
