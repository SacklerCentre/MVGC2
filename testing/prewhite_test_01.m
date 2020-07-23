%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%

n = 10;
p = 20;
w = 0.4;
r = 0.999;
g = 0.3;

m = 40000;

q = 100;
o = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%

A = var_rand(n,p,r,w,'');
V = corr_rand(n,g);
var_info(A,V);

x = varfima_to_tsdata(A,[],[],V,m);
[mosc(1),mosc(2),mosc(3),mosc(4)] = tsdata_to_varmo(x,q,'LWR',[],true,true,1,'');
px = mosc(o);
[Ax,Vx] = tsdata_to_var(x,px,'LWR');
var_info(Ax,Vx);

y = zeros(size(x));
for i = 1:n
	[mosc(1),mosc(2),mosc(3),mosc(4)] = tsdata_to_varmo(x(i,:),q,'LWR',[],false,false,0);
	p = mosc(o);
	[a,v] = tsdata_to_var(x(i,:),p,'LWR');
	fprintf('p = %2d,  rho = %g\n',p,var_specrad(a));
	y(i,:) = genvma(-a,x(i,:));
end

[mosc(1),mosc(2),mosc(3),mosc(4)] = tsdata_to_varmo(y,q,'LWR',[],true,true,1,'');
py = mosc(o);
[Ay,Vy] = tsdata_to_var(y,py,'LWR');
var_info(Ay,Vy);

for i = 1:n
	[mosc(1),mosc(2),mosc(3),mosc(4)] = tsdata_to_varmo(y(i,:),q,'LWR',[],false,false,0);
	p = mosc(o);
	[a,v] = tsdata_to_var(y(i,:),p,'LWR');
	fprintf('p = %2d,  rho = %g\n',p,var_specrad(a));
end
