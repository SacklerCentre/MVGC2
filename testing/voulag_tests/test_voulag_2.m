T     = 10;       % s

tau   = 10;        % ms
theta = 4;        % Hz

dt    = 0.1;      % ms

%%%%%%%%%%%%%%%%%%%%%%%%

tau = tau/1000;
dt  = dt/1000;

K = round(T/dt);
ktau = round(tau/dt);

x = randn(K+1,2);
% x = zeros(K+1,2); x(1:ktau+1,:) = randn(ktau+1,2);
for k = ktau+1:K
	x(k+1,1) = x(k+1,1)+(1-d*dt)*x(k,1) - 2*pi*theta*x(k-ktau,2)*dt;
	x(k+1,2) = x(k+1,2)+(1-d*dt)*x(k,2) + 2*pi*theta*x(k-ktau,1)*dt;
end
t = (0:K)'*dt;

gp_qplot(t,x);
