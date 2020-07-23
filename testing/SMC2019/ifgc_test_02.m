nx = 20;
ny = 15;
m  = 400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
xx = 1:nx;
xy = nx+1:nx+ny;

[V,C0] = corr_rand(n,0.5);

D0 = [C0(xx,xx) zeros(nx,ny)]
D = cell(m,1);
for k = 1:m
	Ck = randn(n)*C0;
	D{k} = Ck(xx,:);
end

DD = zeros((m+1)*nx,(m+1)*n);
for i = 0:m
	for j = 0:i
		if i == j
			DD(i*nx+1:(i+1)*nx,j*n+1:(j+1)*n) = D0;
		else
			DD(i*nx+1:(i+1)*nx,j*n+1:(j+1)*n) = D{i-j};
		end
	end
end

dd = DD(m*nx+1:end,1:m*n);

DD1 = DD(1:m*nx,1:m*n);

D = dd*dd' + D0*D0';

mn = m*n;
S = abs(dd'*inv(D)*dd);

[maxabs(S) 1/m]
