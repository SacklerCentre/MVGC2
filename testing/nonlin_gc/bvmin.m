function [x,y,epsx,epsy] = bvmin(cfb,a,c,xi,eta,kappa,m,mtrans)

if isempty(mtrans), mtrans = m; end

mtot = m+mtrans+1;

u = randn(mtot,2);
epsx = xi*u(:,1);
epsy = kappa*eta*u(:,1) + sqrt(1-kappa*kappa)*eta*u(:,2);

x = epsx;
y = epsy;


z = (x.^cfb(1)).*(y.^cfb(2));

for t = 2:mtot
	x(t) = x(t) + a*x(t-1) + c*z(t-1);
end

m1 = mtrans+2;

x = x(m1:mtot);
y = y(m1:mtot);

if nargout > 2
	epsx = epsx(m1:mtot);
	epsy = epsy(m1:mtot);
end
