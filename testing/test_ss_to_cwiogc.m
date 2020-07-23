n = 6;
m = 20;
rho = 0.9;
g = 1;

%-------------------------------------------------------------

[A,C,K] = iss_rand(n,m,rho);
V = corr_rand(n,g);

ss_info(A,C,K,V,true);

FI = ss_to_cwiogc(A,C,K,V,'in');
FO = ss_to_cwiogc(A,C,K,V,'out');

[FI FO]

GI = zeros(n,1);
GO = zeros(n,1);

for i = 1:n
	r = 1:n; r(i) = [];
	GI(i) = ss_to_mvgc(A,C,K,V,i,r);
end

for i = 1:n
	r = 1:n; r(i) = [];
	GO(i) = ss_to_mvgc(A,C,K,V,r,i);
end

[GI GO]

maxabs(FI-GI)
maxabs(FO-GO)
