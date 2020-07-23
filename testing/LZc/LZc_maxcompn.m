function n = LZc_maxcompn(d,K)

g = 0;
D = 1;
for k = 1:K
	D = D*d;
	g = g+k*D;
end
n = g;
