
nx  = 7;
ny  = 1;
rho = 0.9;
w   = 1.0;
g   = 0.2;
p   = 10;

cmax = 10;

numc = 100;
nums = 100;

%plotm = 'x11';
plotm = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

CON = ones(n); CON(y,x) = 0;

c = linspace(0,cmax,numc)';

gx = zeros(numc,nums);
vx = zeros(numc,nums);
for s = 1:nums
	fprintf('sample %2d of %2d\n',s,nums);

	A = var_rand(CON,p,rho,w,plotm);
	AXY = A(x,y,:);

	V = eye(n);
	V(x,x) = corr_rand(nx,g);

	%gx1 = zeros(numc,1);
	for i = 1:numc
		A(x,y,:) = c(i)*AXY;
		[~,VX] = var2riss(A,V,y,x);
		gx(i,s) = sum(log(diag(VX))) - logdet(VX);
		vx(i,s) = logdet(VX);
		%[AA,CC,KK] = var_to_ss(A,V);
		%[~,VX1] = ss2iss(AA,CC(x,:),KK*V*KK',V(x,x),KK*V(:,x));
		%gx1(i) = sum(log(diag(VX1))) - logdet(VX1);
	end
end

gp_qplot(c,gx,'','unset key\nset yr[0:*]\nset ylab "-log|R|" rot\nset xlab "strength of latent influence"\nset title "Residuals (log-)generalised correlation"\nset grid','epsl',[],14,'~/tmp/gencorr');
gp_qplot(c,vx,'','unset key\nset yr[*:*]\nset ylab "log|V|" rot\nset xlab "strength of latent influence"\nset title "Residuals (log-)generalised variance"\nset grid','epsl',[],14,'~/tmp/genvar');
