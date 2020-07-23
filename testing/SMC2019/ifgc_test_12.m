a    = 0.95;
b    = -0.9;
c    = 2;
t    = 20;
g    = 0.5;
kmax = 50;

tol = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx  = 1;
ny  = 1;

n = nx+ny;
x = 1:nx;
y = nx+1:n;

VARA(:,:,1) = [a 0; 0 b];
VARA(:,:,t) = [0 c; 0 0];

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
