T     = 2;       % s

dect  = 100;     % ms

theta = 60;      % Hz

dt    = 0.1;      % ms

%%%%%%%%%%%%%%%%%%%%%%%%

NODE = [dect theta]
CON = zeros(0,4);
V = 1;

[X,t] = voulag(NODE,CON,V,dt,T,0,false);

gp_qplot(t',X');

fs = 1/(dt/1000);
window = 10000;
[S,f,nwobs,noobs,nwins] = tsdata_to_cpsd(X,fs,window	,[],[],true);

nwobs
noobs
nwins

gpcmds = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead\nset logs x',theta,theta);
gp_qplot(f,log(S),[],gpcmds);
