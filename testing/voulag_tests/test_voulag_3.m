T     = 10;       % s

dect  = 100;     % ms

theta = 60;      % Hz

dt    = 0.1;      % ms

%%%%%%%%%%%%%%%%%%%%%%%%

dt  = dt/1000;
dect = dect/1000;

fs = 1/dt;

K = round(T/dt);

d = 1-dt/dect;
o = 2*pi*dt*theta;

z = randn(K+1,2)*sqrt(dt);
x = zeros(K+1,2);
for k = 1:K
	x(k+1,1) = z(k+1,1) + d*x(k,1) - o*x(k,2);
	x(k+1,2) = z(k+1,2) + d*x(k,2) + o*x(k,1);
end
t = (0:K)'*dt;
gp_qplot(t,x);

% A = [1-d*dt -2*pi*theta*dt; 2*pi*theta*dt 1-d*dt]
% specnorm(A)

window = 2000;
[S,f,nwobs,noobs,nwins] = tsdata_to_cpsd(x',fs,window	,[],[],true);

nwobs
noobs
nwins

gpcmds = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead\nset logs x',theta,theta);
gp_qplot(f,log(S),[],gpcmds);
