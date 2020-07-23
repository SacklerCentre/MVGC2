T     = 20;       % s

fs    = 1000;     % Hz

theta = 40;      % Hz
stheta = 0.0;

sphi = 1;

%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/fs;
m  = round(T*fs);
T  = m/fs;
t  = (0:m)'/fs;

omega = theta*(1 + stheta*dt*randn(m+1,1));

phi = sphi*dt*randn(m+1,1);
for i = 1:m
	phi(i+1) = phi(i) + phi(i+1);
end

x = sin(2*pi*omega.*(t-phi));
gp_qplot(t,x);

window = 1000;
[S,f,nwobs,noobs,nwins] = tsdata_to_cpsd(x',fs,window,[],[],true);

nwobs
noobs
nwins

gpcmds = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead',theta,theta);
gpcmds = [gpcmds '\nset logs x'];
gp_qplot(f,log(S),[],gpcmds);
