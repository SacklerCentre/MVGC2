function [gama,iters] = corr_rand_gama(n,gtarget,S,tol,imax,verb)

	if nargin < 4 || isempty(tol),  tol  = 1e-08; end
	if nargin < 5 || isempty(imax), imax = 100;   end
	if nargin < 6 || isempty(verb), verb = false; end

	g = zeros(S,1);

	pmax = 2;
	gmean = 0;
	k = 1;
	while gmean < gtarget
		pmax = pmax/2;
		gmean = meang(n,pmax,S);
		if verb, fprintf('init %2d : pmax = %7.4f, g = %8.4f\n',k,pmax,gmean); end
		k = k+1;
	end
	if verb, fprintf('\n'); end

	pmin = 0;
	gama = Inf;
	dgama = Inf;
	iters = 1;
	while iters <= imax && abs(dgama) > tol
		ogama = gama;
		gama = (pmin+pmax)/2;
		dgama = gama-ogama;
		gmean = meang(n,gama,S);
		if verb, fprintf('trial %3d : gama = %7.4f, g = %8.4f, dgama = % g\n',iters,gama,gmean,dgama); end
		if gmean > gtarget
			pmin = gama;
		else
			pmax = gama;
		end
		iters = iters+1;
	end

end

function gmean = meang(n,gama,S)

	for s = 1:S
		[Q,R] = qr(randn(n));
		M = Q*diag(sign(diag(R))); % M orthogonal
		v = gamrnd(gama,2,n,1);
		g(s) = sum(log(diag(M*diag(v)*M'))-log(v));
	end
	gmean = mean(g);
end
