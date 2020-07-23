nx   = 4;
ny   = 5;
nz   = 0;
rho  = 0.9;
P    = 7;
p    = 20;
w    = 1.0;
g    = 0.5;

kmax = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
u = randperm(n);
x = u(1:nx)
y = u(nx+1:nx+ny)
z = u(nx+ny+1:n)

VARA = zeros(n,n,p);
VARA(:,:,1:P) = var_rand(n,P,rho,w);

if kill2
	for k = 2:p
		A1Ak1 = VARA(:,:,1)*VARA(:,:,k-1);
		VARA(x,y,k) = -A1Ak1(x,y);
	end
else % kill1
	VARA(x,y,:) = 0;
end

specnorm(VARA)

V = corr_rand(n,g);

[A,C,K,info] = var_to_ss(VARA,V);

F = ss_to_fhgc(A,C,K,V,x,y,kmax);

F1 = ss_to_mvgc(A,C,K,V,x,y)

[F(1) F(2)]

k = (1:kmax)';

gp_qplot(k,F,'multi-step GC');
