nx = 4;
ny = 5;
nz = 2;
m = 11;
rho = 0.9;
g = 1;

%-------------------------------------------------------------

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;

[A,C,K] = iss_rand(n,m,rho,1);
V = corr_rand(n,g);

ss_info(A,C,K,V,true);

[F,P] = ss_to_mvgc(A,C,K,V,x,y)

return

P11 = P(1:n,1:n)
P22 = P(n+1:2*n,n+1:2*n)
P12 = P(1:n,n+1:2*n)
