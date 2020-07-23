nx   = 4;
ny   = 5;
nz   = 3;
p    = 20;
rho  = 0.9;
w    = 1.0;
g    = 0.5;
kmax = 50;

t = 5:8;
tol = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
u = randperm(n);
x = u(1:nx)
y = u(nx+1:nx+ny)
z = u(nx+ny+1:n)

VARA = var_rand(n,p,rho,w);
tt = 1:p;
tt(t) = [];
VARA(:,:,tt) = 0;
VARA = specnorm(VARA,rho);

V = corr_rand(n,g);

[A,C,K,ssinfo] = var_to_ss(VARA,V);

[F0,~,dF0] = ss_to_ffgc(A,C,K,V,x,y,kmax);
[G2,~,dF2] = ss_to_msgc(A,C,K,V,x,y,kmax);

[dF0 dF2]

G0 = [NaN;diff(F0)];
F2 = cumsum(G2);

k = (1:kmax)';
gp_qplot(k,[F0 F2],{'FF','MS'},'set title "CUMULATIVE GC"\nset key top left');
gp_qplot(k,[G0 G2],{'FF','MS'},'set title "PER-STEP GC"\nset key top right');
