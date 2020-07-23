%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usemex = false;

simt  = 10;   % total simulation time (secs)
sett  = 0;    % setle time, to allow for initial transient to decay (secs)

dt = 0.1;     % simulation time step (ms)
ds = 1;       % downsample factor (so sample frequency after downsampling is 1000/dt/ds

n  = 9; % number of nodes
G = [   % causal graph specification
	1 3;
	1 4;
	2 4
	3 4;
	4 2;
	5 3
	6 5;
	7 6;
	8 9;
	9 7
];

mdect = 100; % node decay: mean (ms)
sdect =  10; % node decay: std  (ms)
mclag =  20; % causal lag: mean (ms)
sclag =   5; % causal lag: std  (ms)
mcfbs =   1; % causal feedback: mean (1/secs)
scfbs =   0.1; % causal feedback: std  (1/secs)

pncfb = 0.25; % probability of a negative causal feedback

oscs  = [   % oscillator nodes, and frequencies (Hz)
	4 60;
	6 20
];

mseed = 0; % model random seed (0 for unseeded)
xseed = 0; % time-series random seed (0 for unseeded)

gpterm = 'x11';

%%%%%%%%%%%%%%%%%% Set up causal model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = size(G,1); % number of connections

FS = 1000/dt; % simulation frequency (Hz)

rstate = rng_seed(mseed);
dect = gamrnd((mdect^2)/(sdect^2),(sdect^2)/mdect,n,1); % Gamma-distributed node decay
clag = gamrnd((mclag^2)/(sclag^2),(sclag^2)/mclag,c,1); % Gamma-distributed causal lag
cfbs = gamrnd((mcfbs^2)/(scfbs^2),(scfbs^2)/mcfbs,c,1); % Gamma-distributed causal feedback
cfbs = (2*(rand(c,1)>pncfb)-1).*cfbs;
V  = corr_rand(n,1); % residuals covariance matrix; so dW = V*sqrt(dt)
rng_restore(rstate);

NODE = [dect zeros(n,1)];
NODE(oscs(:,1),2) = oscs(:,2)

CON = [G clag cfbs]  % the causal connection configuration

%%%%%%%%%%%%%%%%%% Generate lagged OU data and downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rstate = rng_seed(xseed);
%tic
[X,t] = voulag(NODE,CON,V,dt,simt,sett,usemex);
%toc
rng_restore(rstate);

[Y,fs] = downsample(X,ds,FS);
T      = downsample(t,ds,FS);

window = 1000;
[S,f] = tsdata_to_cpsd(Y,fs,window,[],[],true);

%%%%%%%%%%%%%%%%%% Plot time-series and spectral power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(gpterm)
	figure(1);clf;
	sgtitle({['sample frequency = ' num2str(fs) 'Hz'],''});

	subplot(2,1,1);
	plot(T',Y');
	xlabel('time (seconds)');
	ylabel('signal');

	subplot(2,1,2);
	plot(f,log(S));
	xlabel('frequency (Hz)');
	ylabel('power (dB)');
	xlim([f(1) f(end)]);
	set(gca,'XScale','log');
else
	gp_qplot(T',Y',[],['set title "sample frequency = ' num2str(fs) 'Hz"\nset xlabel "time (seconds)"\nset ylabel "signal" rot'],gpterm);
	gp_qplot(f,log(S),[],['	set title "sample frequency = ' num2str(fs) 'Hz"\nset xr [1:*]\n\nset xlabel "frequency (Hz)"\nset ylabel "power (dB)" rot\nset logs x\nset grid'],gpterm);
end
