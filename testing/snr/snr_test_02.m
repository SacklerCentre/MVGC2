%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx        = 1;
ny        = 9;
nz        = 0;

p         = 10;     % AR model order
rho       = 0.9;    % AR spectral radius
w         = 0.5;    % var coefficients decay weighting factor
d         = 0.75;   % fractional integration exponent
r         = 20;     % fractional integration lags
rmi       = 0.1;    % residuals log-generalised correlation
nrmi      = 2.0;    % noise log-generalised correlation (multi-information)

lmin = 0;
lmax = 5;
lres = 100;

fres = 600;
fs   = 200;

if ~exist('seed', 'var'), seed = 0; end % random seed (0 for unseeded)

gpterm = 'epsl';
gpname = fullfile(tempdir,'snrtest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed);

f = sfreqs(fres,fs);

n = nx+ny+nz;
prm = randperm(n);
x = prm(1:nx);
y = prm(nx+1:nx+ny);
z = prm(nx+ny+1:n);

while true
	VARA1 = var_rand(n,p,rho,w);

	c = fracint_coeffs(d,r,true);

	VARA  = mvconv(cat(3,eye(n),-VARA1),c);
	VARA = -VARA(:,:,2:end);

	% VARAref = arfilt(VARA1,c);
	% maxabs(VARAref-VARA)
	% specnorm(VARA)

	V = corr_rand(n,rmi);

	[A,C,K] = var_to_ss(VARA,V);

	LV = chol(V,'lower');
	KL = K*LV;
	KV = K*V;
	Q = KL*KL';

	O = dlyap(A,Q);
	LO = chol(O,'lower');
	CLO = C*LO;
	W = CLO*CLO' + V; % cov matrix of reference signal

	N = corr_rand(n,nrmi); % cov matrix of noise signal
	N = mean(diag(W)./diag(N))*N;

	I = cov_to_mvmi(W,x,y); % (x,y|z) partial correlation for ref signal
	J = cov_to_mvmi(N,x,y); % (x,y|z) partial correlation for noise signal

	[I J]
	if J > I, break; end
end

ss_info(A,C,K,V);

S = ss_to_cpsd(A,C,K,V,fres);
Sa = zeros(fres+1,n);
for i = 1:n, Sa(:,i) = 20*log10(squeeze(S(i,i,:))); end
gp_qplot(f,Sa,[],'set logs x\nunset key\nset xlabel "frequency (Hz)\nset ylabel "dB"\n',gpterm,0.8,[],[gpname '_PSD']);

res = input('\nContinue? (y/N) ','s');
if isempty(res) || ~strcmpi(res,'y'), return; end

F = ss_to_mvgc(A,C,K,V,x,y);

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

	CCLO = CC*LO;
	WW = CCLO*CCLO' + VV; % cov matrix of composite signal

	II(i) = cov_to_mvmi(WW,x,y);
end
fprintf('\n\n');

%SNR = lams2;
SNR = 20*log10(lams); % dB

hline1 = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead lc "red"\n',F,F);
htext1 = sprintf('set label "pure signal GC" at graph 0.1,first %g left front boxed\n',F);
hline2 = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead lc "red"\n',I,I);
htext2 = sprintf('set label "pure signal MI" at graph 0.1,first %g left front boxed\n',I);
hline3 = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead lc "green"\n',J,J);
htext3 = sprintf('set label "pure noise MI" at graph 0.1,first %g left front boxed\n',J);
cmds = 'unset key\nset xlabel "high noise <--- SNR (dB)---> low noise"\nunset ylabel\nset style textbox opaque noborder\n';

yr1 = sprintf('set yr [0:%g]\n',1.05*F);
yr2 = sprintf('set yr [%g:%g]\n',0.9*min(min(II),min(I,J)),1.05*max(max(II),max(I,J)));

gp_qplot(SNR,FF,[],['set title "Granger causality"\n'  cmds yr1 hline1 htext1],gpterm,0.8,[],[gpname '_GC']);
gp_qplot(SNR,II,[],['set title "Mutual information"\n' cmds yr2 hline2 htext2 hline3 htext3],gpterm,0.8,[],[gpname '_MI']);

%{
function B = arfilt(A1,c)

	[n,~,p] = size(A1);
	r = length(c)-1;
	q = p+r;

	A = cat(3,eye(n),-A1);

	for k = 0:q
		i1 = max(0,k-r);
		i2 = min(p,k);
		Bk = zeros(n);
		for i = i1:i2
			Bk = Bk + A(:,:,i+1)*c(k-i+1);
		end
		B(:,:,k+1) = Bk;
	end

	B = -B(:,:,2:end);

end
%}
