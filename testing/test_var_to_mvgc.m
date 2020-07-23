nx = 4;
ny = 5;
nz = 2;
p = 3;
rho = 0.9;
w = 0.5;
g = 1;

rng_seed(21312);

%-------------------------------------------------------------

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;

VARA = var_rand(n,p,rho,w);
V = corr_rand(n,g);

[A,C,K,info] = var_to_ss(VARA,V);

[F,P] = ss_to_mvgc(A,C,K,V,x,y)

P11 = P(1:n,1:n)
P22 = P(n+1:2*n,n+1:2*n)
P12 = P(1:n,n+1:2*n)
