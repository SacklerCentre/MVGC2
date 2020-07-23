r   = 20;
nx  = 5;
ny  = 14;
rho = 0.95;
g   = 1;
m   = 100;

tol = sqrt(eps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

[A,C,K] = iss_rand(n,r,rho);

[V,C0] = corr_rand(n,g);

D0 = C0(1:nx,:);
D = cell(m,1);
Ak1 = eye(r);
for k = 1:m
	Ck = C*Ak1*K*C0; % Bk = C*A^{k-1}*K is the k-th MA coefficients matrix
	D{k} = Ck(x,:);
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

KVL = K*chol(V,'lower');
[~,VR,rep] = ss2iss(A,C(x,:),KVL*KVL',V(x,x),K*V(:,x));
if sserror(rep), return; end
LDVR = logdet(VR);

F0 = LDVR-logdet(V(x,x))
F = nan(m,1);
for k = 1:m
	DDk = DD(1:(k+1)*nx,1:(k+1)*n);
	F(k) = (k+1)*LDVR-logdet(DDk*DDk');
	if k == 1
		FDk = F(k)-F0;
	else
		FDk = F(k)-F(k-1);
	end
	fprintf('%4d  %g\n',k,FDk);
	if FDk < tol
		break;
	end
end
mm = (1:m+1)';

[F0f,Ff] = ifgc(A,C,K,V,x,y,m,tol);

fprintf('\nitterations = %d, diff = %g\n\n',k,maxabs([F0;F]-[F0f;Ff]));

gp_qplot(mm,[F0;F],{'GC'},'set key bottom right');
