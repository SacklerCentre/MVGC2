%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usemex = false;
mseed  = 0; % model random seed (0 for unseeded)
xseed  = 0; % time-series random seed (0 for unseeded)

gpterm = '';

tsim  = 80;  % total simulation time (secs)

fs = 10000;  % simulation frequency (Hz)
ds = 20;      % downsample factor (so sample frequency after downsampling is fs/ds

%  to  from
G = [        % causal graph specification
	1  3;
	1  4;
	2  4;
	3  4;
	4  2;
	5  3;
	6  5;
	7  6;
	8  9;
	9  7
];
n = max(G(:,1)); % number of nodes
c = size(G,1);    % number of connections

on = 2;
oc = on*n;
OG = zeros(oc,2);
for k = 1:on;
	OG((k-1)*n+1:k*n,:) = [(1:n)' k*ones(n,1)];
end

fdecm  =   4;      % node decay frequency: mean (Hz)
fdecs  =   1;      % node decay frequency: std  (Hz)

tlagm  =  20/1000; % causal lag: mean (ms)
tlags  =   5/1000; % causal lag: std  (ms)
fcaum  =   3.5;    % causal feedback: mean (Hz)
fcaus  =   0.5;    % causal feedback: std  (Hz)
pncau  =   0.25;   % probability of a negative causal feedback

fodca  =   0.2;    % alpha node decay frequency
fosca  =  12;      % alpha node oscillator frequency
fodcd  =   0.1;    % delta node decay frequency
foscd  =   3;      % delta node oscillator frequency

focama =   0.5;    % alpha node causal feedback mean
focasa =   0.02;   % alpha node causal feedback std

focamd =   0.3;    % delta node causal feedback mean
focasd =   0.02;   % delta node causal feedback std

%%%%%%%%%%%%%%%%%% Set up causal model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rstate = rng_seed(mseed);
fdec  = gamrnd((fdecm^2)/(fdecs^2),(fdecs^2)/fdecm,n,1);     % Gamma-distributed node decay
fcau  = gamrnd((fcaum^2)/(fcaus^2),(fcaus^2)/fcaum,c,1);     % Gamma-distributed causal feedback
fcau  = (2*(rand(c,1)>pncau)-1).*fcau;                       % excitatory/inhibitory
tlag  = gamrnd((tlagm^2)/(tlags^2),(tlags^2)/tlagm,c,1);     % Gamma-distributed causal lag
focaa = gamrnd((focama^2)/(focasa^2),(focasa^2)/focama,n,1); % alpha node Gamma-distributed causal feedback
focad = gamrnd((focamd^2)/(focasd^2),(focasd^2)/focamd,n,1); % delta node Gamma-distributed causal feedback
V     = corr_rand(n,1);                                      % residuals covariance matrix (left-Cholesky factor
rng_restore(rstate);

foca = [focaa;focad];

NODE  = fdec;
CONX  = [G fcau tlag];                     % the causal connection configuration
ONODE = [fodca fosca; fodcd foscd];
OCONX = [OG foca];

%%%%%%%%%%%%%%%%%% Generate lagged OU data and downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rstate = rng_seed(xseed);
tic
[X,t] = voulag(tsim,fs,NODE,CONX,ONODE,OCONX,V,usemex);
toc
rng_restore(rstate);

[tX,fs] = downsample([t;X],ds,fs);
t = tX(1,:);
X = tX(2:end,:);

%[S,f] = periodogram(X',[],[],fs);
window = 2^12;
fres   = 2^12;
[S,f]  = tsdata_to_cpsd(X,fs,window,[],fres,true);

SdB = 20*log10(S);

%%%%%%%%%%%%%%%%%% Plot time-series and spectral power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(gpterm)
	figure(1);clf;
	sgtitle({['sample frequency = ' num2str(fs) 'Hz'],''});

	subplot(2,1,1);
	plot(t',X');
	xlabel('time (seconds)');
	ylabel('signal');

	subplot(2,1,2);
	plot(f,SdB);
	xlabel('frequency (Hz)');
	ylabel('power (dB)');
	xlim([1 f(end)]);
	set(gca,'XScale','log');
	xline(fosca);
	xline(foscd);
	grid on;
else
	gp_qplot(t',X',[],[],gpterm);
	gpcmds = '';
	gpcmds = [gpcmds sprintf('\nset arrow from first %g,graph 0 to first %g,graph 1 nohead ls 1 lw 2',fosca,fosca)];
	gpcmds = [gpcmds sprintf('\nset arrow from first %g,graph 0 to first %g,graph 1 nohead ls 2 lw 2',foscd,foscd)];
	gpcmds = [gpcmds '\nset xr [1:500]'];
	gpcmds = [gpcmds '\nset logs x'];
	gp_qplot(f,SdB,[],gpcmds,gpterm);
end

%%%%%%%%%%%%%%%%%% Granger causalities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regmode = 'LWR';
varmomax = 30;
ssmomax  = 30;
alpha = 0.05;
mhtc = 'BONFERRONI';

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,varmomax,'LWR',[],[],2);

varmo = molrt;

[A,V] = tsdata_to_var(X,varmo,regmode);
info = var_info(A,V);
assert(~info.error,'VAR error(s) found - bailing out');

[F1,pval] = var_to_pwcgc(A,V,X,regmode);
assert(~isbad(F1,false),'GC estimation failed');

sig = significance(pval.FT,alpha,mhtc);

GG = zeros(n);
for k = 1:c, GG(G(k,1),G(k,2)) = 1; end

pf = 2*varmo; % Bauer recommends 2 x VAR model order

[ssmo,ssmomax] = tsdata_to_sssvc(X,pf,[],3);
assert(ssmo > 0,'selected zero model order! GCs will all be zero!');
if ssmo >= ssmomax, fprintf(2,'*** WARNING: selected SS maximum model order (may have been set too low)\n'); end

[A,C,K,V] = tsdata_to_ss(X,pf,ssmo);
info = ss_info(A,C,K,V);
assert(~info.error,'SS error(s) found - bailing out');

F2 = ss_to_pwcgc(A,C,K,V);
assert(~isbad(F2,false),'GC estimation failed');

maxF = 1.1*nanmax([F1(:);F2(:)]);
pdata = {GG sig;F1 F2};
ptitle = {'ground truth' 'F-test (VAR)'; 'PWCGC (VAR)' 'PWCGC (SS)'};
maxp = [1 1; maxF maxF];
plot_gc(pdata,ptitle,[],maxp,4);
