r = 2;
nx = 1;
ny = 1;
a = 0.7;
b = 0.8;
c = 1;
rho = 0.3;

kmax   = 100;

%tol = eps;
tol = sqrt(eps);

S = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


v = linspace(v1,v2,S)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

F = nan(S,1);

for i = 1:S
	eval(sprintf('%s = v(%d);',vname,i));

	VARA = [a c; 0 b];
	V = [1 rho; rho 1];
	[A,C,K] = var_to_ss(VARA,V);

	[FF,k,dF,] = ifgc(A,C,K,V,x,y,kmax,tol);
	if dF < tol
		F(i) = FF(end);
	end
end

gpcmds = sprintf('unset key\nset xlabel "variable: %s"',vname);
gp_qplot((1:S)',F,[],gpcmds);
