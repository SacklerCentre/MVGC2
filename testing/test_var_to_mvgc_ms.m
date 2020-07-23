nx  = 3;
ny  = 6;
nz  = 5;
p   = 10;
rho = 0.95;
w   = 1.0;
g   = 0.5;

q = 40;

%rng_seed(21312);

%-------------------------------------------------------------

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;

VARA = var_rand(n,p,rho,w);
%VARA = zeros(n,n,p); VARA(:,:,p) = randn(n); VARA = var_specrad(VARA,rho);
V    = corr_rand(n,g);

F1 = var_to_mvgc(VARA,V,x,y)
F = var_to_mvgc_ms(VARA,V,x,y,q)

gp_qplot((1:q)',F)
