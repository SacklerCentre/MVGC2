r = 2;
nx = 1;
ny = 1;
a = 0.7;
b = 0.8;
c = 5;
rho = 0.3;

t = 10; % lag

kmax = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

VARA = zeros(2,2,t);
VARA(:,:,1) = [a 0; 0 b];
VARA(:,:,t) = [0 c; 0 0];
V = [1 rho; rho 1];
[A,C,K] = var_to_ss(VARA,V);

FH = fhgc(A,C,K,V,x,y,kmax);
FI = ifgc(A,C,K,V,x,y,kmax);

if perstep
	FI = [FI(1);diff(FI)];
	gpcmds = 'set key top right Left rev\nset title "Per-step GC"';
else
	FH = cumsum(FH);
	gpcmds = 'set key bottom right Left rev\nset title "Cumulative GC"';
end

k = (1:kmax)';

gpcmds = sprintf('%s\nset xr [1:%d]\nset arrow from first %d,graph 0 to first %d,graph 1 nohead',gpcmds,kmax,t,t);
gp_qplot(k,[FH FI],{'finite horizon','total'},gpcmds);
