n = 7;
m = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,C0] = corr_rand(n,1);

for k = 1:m
	C{k} = randn(n)*C0;
end

CC = zeros((m+1)*n,(m+1)*n);
for i = 0:m
	ii = i*n+1:(i+1)*n;
	for j = 0:i
		if i == j
			CC(ii,ii) = C0;
		else
			jj = j*n+1:(j+1)*n;
			CC(ii,jj) = C{i-j};
		end
	end
end

EE = CC*CC';

log(det(EE))/(m+1)
2*(log(det(CC))/(m+1))
2*(log(prod(diag(CC)))/(m+1))
2*log(det(C0))
log(det(V))

return


cc = CC(m*n+1:(m+1)*n,1:m*n);

%DD1 = CC(1:m*n,1:m*n);

Delta = cc*cc' + C0*C0';

S = cc'*inv(Delta)*cc;

maxabs(S)
