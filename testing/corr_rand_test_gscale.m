
ens  = [20 40 60 80 100 120 140 160 180 200]';

pees = linspace(2,20,19)';

S = 1000;
P = length(pees);
N = length(ens);

a = zeros(N,1);
b = zeros(N,1);

for i = 1:N
	n = ens(i);
	fprintf('trial %2d of %2d : n = %g\n',i,N,n);
	g = zeros(S,P);
	for k = 1:P
		p = pees(k);
		fprintf('\tsample %2d of %2d : p = %g',k,P,p);
		for s = 1:S
			[Q,R] = qr(randn(n));
			M = Q*diag(sign(diag(R))); % M orthogonal
			v = realpow(abs(randn(n,1)),p);
			g(s,k) = sum(log(diag(M*diag(v)*M'))-log(v));
		end
		fprintf('\n')
	end
	gmean = mean(g)';
	pp = [pees ones(P,1)];
	c = pp\gmean;
	a(i) = c(1);
	b(i) = c(2);
end

gp_qplot(ens,a);
