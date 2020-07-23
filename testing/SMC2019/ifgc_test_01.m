nx = 2;
ny = 3;
m  = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
xx = 1:nx;
xy = nx+1:nx+ny;

[V,C0] = corr_rand(n,0.5);

E0 = C0(xx,xx);
F0 = zeros(nx,ny);
for k = 1:m
	Bk = randn(n);
	Ck = Bk*C0;
	E{k} = Ck(xx,xx);
	F{k} = Ck(xx,xy);
end

ee1 = zeros(m*nx,nx);
ff1 = zeros(m*nx,ny);
for k = 1:m
	ee1((k-1)*nx+1:k*nx,:) = E{k};
	ff1((k-1)*nx+1:k*nx,:) = F{k};
end

EE = zeros((m+1)*nx,(m+1)*nx);
for i = 0:m
	for j = 0:i
		if i == j
			EE(i*nx+1:(i+1)*nx,j*nx+1:(j+1)*nx) = E0;
		else
			EE(i*nx+1:(i+1)*nx,j*nx+1:(j+1)*nx) = E{i-j};
		end
	end
end

FF = zeros((m+1)*nx,(m+1)*ny);
for i = 0:m
	for j = 0:i
		if i == j
			FF(i*nx+1:(i+1)*nx,j*ny+1:(j+1)*ny) = F0;
		else
			FF(i*nx+1:(i+1)*nx,j*ny+1:(j+1)*ny) = F{i-j};
		end
	end
end

ee = EE(nx+1:end,1:nx);
ff = FF(nx+1:end,1:ny);

eeb = EE(m*nx+1:end,1:m*nx);
ffb = FF(m*nx+1:end,1:m*ny);

EE1  = EE(1:m*nx,1:m*nx);
FF1  = FF(1:m*nx,1:m*ny);
ee1  = EE1(nx+1:end,1:nx);
ff1  = FF1(nx+1:end,1:ny);
eeb1 = EE1((m-1)*nx+1:end,1:(m-1)*nx);
ffb1 = FF1((m-1)*nx+1:end,1:(m-1)*ny);

EE2  = EE(1:(m-1)*nx,1:(m-1)*nx);
FF2  = FF(1:(m-1)*nx,1:(m-1)*ny);
ee2  = EE2(nx+1:end,1:nx);
ff2  = FF2(nx+1:end,1:ny);
eeb2 = EE1((m-2)*nx+1:end,1:(m-2)*nx);
ffb2 = FF1((m-2)*nx+1:end,1:(m-2)*ny);

G = EE1*EE1' + FF1*FF1' + ff*ff'

G1 = EE2*EE2' + FF2*FF2' + ff1*ff1'

D = E0*E0' + eeb1*eeb1' + ffb*ffb'

g = EE2*eeb1' + FF2*ffb1' + ff1*F{m}'

H = g*inv(D)*g'

G1 - H

M = g'*inv(G1)*g

D-M
