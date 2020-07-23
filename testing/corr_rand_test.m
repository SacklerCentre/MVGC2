gpterm = 'png';

p = linspace(.1,20,100)';
pfit = 4;

S = 1000;
P =length(p);

g = zeros(S,P);
for k = 1:P
	fprintf('sample %2d of %2d : p = %g',k,P,p(k));
	for s = 1:S
		[Q,R] = qr(randn(n));
		M = Q*diag(sign(diag(R))); % M orthogonal
		v = realpow(abs(randn(n,1)),p(k));
		g(s,k) = sum(log(diag(M*diag(v)*M'))-log(v));
	end
	fprintf('\n')
end

gmean = mean(g)';

y = gmean(p >= pfit);
x = p(p >= pfit);
m = length(x);


xx = [x ones(m,1)];
c = xx\y;
a = c(1);
b = -c(2);

gfit1 = a*p-b;
gfit2 = a*(p-2);
gfit3 = a*p;

b/a

gunif = multiinfo(n,true)

leg = {sprintf('mean(g) : n = %d',n), sprintf('linear fit: a = %g, b = %g',a,b)};
guline = sprintf('set arrow from graph 0,first %g to graph 1,first %g nohead\n',gunif,gunif);
cmds = ['set grid\nset xzeroaxis lt 1\nset xlabel "exponent"\nset ylabel "mean multi-info" rot\nset key top left Left rev\n' guline];
gp_qplot(p,[gmean gfit1 gfit2 gfit3],leg,cmds,gpterm);
