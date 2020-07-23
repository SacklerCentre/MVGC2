r   = 11;
nx  = 4;
ny  = 5;
rho = 0.95;
g   = 1;
kmax   = 100;

tol = [];

%rng_seed(867817212);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

[A,C,K] = iss_rand(n,r,rho);

V = corr_rand(n,g);

ssinfo = ss_info(A,C,K,V);

[F1,~,dF1] = ss_to_ffgc(A,C,K,V,x,y,kmax,tol);
[G2,~,dF2] = ss_to_msgc(A,C,K,V,x,y,kmax,tol);

[dF1 dF2]

G1 = [NaN;diff(F1)];
F2 = cumsum(G2);

k = (1:kmax)';
gp_qplot(k,[F1 F2],{'FF','MS'},'set title "CUMULATIVE GC"\nset key top left');
gp_qplot(k,[G1 G2],{'FF','MS'},'set title "PER-STEP GC"\nset key top right');
