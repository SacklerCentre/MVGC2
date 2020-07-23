r   = 7;
nx  = 5;
ny  = 4;
rho = 0.9;
g   = 1;
m   = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;

[A,C,K] = iss_rand(n,r,rho);

[V,C0] = corr_rand(n,g);

D0 = C0(1:nx,:);
D = cell(m,1);
Ak1 = eye(r);
for k = 1:m
	Ck = C*Ak1*K*C0; % Bk = C*A^{k-1}*K is the k-th MA coefficients matrix
	D{k} = Ck(1:nx,:);
	Ak1 = A*Ak1;
end

DD = zeros((m+1)*nx,(m+1)*n);
for i = 0:m
	ii = i*nx+1:(i+1)*nx;
	for j = 0:i-1
		jj = j*n+1:(j+1)*n;
		DD(ii,jj) = D{i-j};
	end
	jj = i*n+1:(i+1)*n;
	DD(ii,jj) = D0;
end

ldk = zeros(m+1,1);
ldk(1) = NaN;
adk = zeros(m+1,1);
adk(1) = NaN;
for k = 1:m
	DDk = DD(1:(k+1)*nx,1:(k+1)*n);
	ldk(k+1) = logdet(DDk*DDk')/k;

	dk = DD(k*nx+1:(k+1)*nx,1:k*n);
	DELk = dk*dk' + D0*D0';
	adk(k+1) = logdet(DELk);
end
mmm = (0:m)';

gp_qplot(mmm,[ldk adk],{'actual','approx'},'set key bottom right');
