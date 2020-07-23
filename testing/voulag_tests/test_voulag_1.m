mseed = 120065;
xseed = 981091;

% ds = 100;

rstate = rng_seed(mseed);
dect  = [30 40 50 35 10 21 25 15 20]';
n     = length(dect);
dt    = 0.01;
s     = 2;
CON   = [ % tnet9 !
	1 3 20 20*s;
	1 4 33  8*s;
	2 4 30 10*s
	3 4 20 20*s;
	4 2 33  8*s;
	5 3 30 10*s
	6 5 20 20*s;
	7 6 33  8*s;
	8 9 30 10*s
	9 7 30 10*s
];

V     = corr_rand(n,1);
rng_restore(rstate);

simt  = 10;
sett  = 0;

rstate = rng_seed(xseed);
tic
[X,t,Z] = voulag(CON,dect,V,dt,simt,sett,true);
toc
rng_restore(rstate);

%{
rstate = rng_seed(xseed);
tic
[Y,t,W] = voulag(CON,dect,V,dt,simt,sett,false);
toc
rng_restore(rstate);
fprintf('\nmaxabs     = %.30f\n',  maxabs(X-Y));
%}

fprintf(  'noise sdev = %.30f\n',  mean(std(Z')));
fprintf(  'proc  sdev = %.30f\n\n',mean(std(X')));

FS = 1000/dt;
[Y,fs] = downsample(X,ds,FS);
T      = downsample(t,ds,FS);

gp_qplot(T',Y',[],['set title "sample frequency = ' num2str(fs) 'Hz"\nset xlabel "time (seconds)"\nset ylabel "signal" rot']);

[S,f,nwobs,noobs,nwins] = tsdata_to_cpsd(Y,fs,[],[],[],true);
gp_qplot(f,log(S),[],['set title "sample frequency = ' num2str(fs) 'Hz"\nset xlabel "frequency (Hz)"\nset ylabel "power (dB)" rot\nset logs x']);
