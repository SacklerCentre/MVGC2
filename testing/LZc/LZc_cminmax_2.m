function [cmin,cmax] = LZc_cminmax_2(d,nn) % max complexity for strings of length n, alphabet size d

cmin = floor((sqrt(8*sum(nn)+1)-1)/2);

g  = 0; % g_k(d)     = d + 2*d^2 + ... + k*d^k
f  = 0; % f_{k-1}(d) = d + d^2 + ... + d^{k-1}
dk = 1; % d^k

n = nn(1);
k = 1;
while true
	dk = dk*d;
	n0 = g; % the D = k' = 0 string length
	g = g+k*dk;
	if g > n % overshot: decrement k
		break
	end
	f = f+dk;
	k = k+1;
end
%cmax = f + floor((n-n0)/k)-1; % f_{k-1}(d) + D, D = floor((n-g_{k-1}(d))/k)
k = k-1;

n = nn(1)+nn(2);
K = 1;
while true
	dk = dk*d;
	n0 = g; % the D = K' = 0 string length
	g = g+K*dk;
	if g > n % overshot: decrement K
		break
	end
	f = f+dk;
	K = K+1;
end
cmax = f + floor((n-n0-k)/K); % f_{K-1}(d) + D, D = floor((n-g_{K-1}(d))/K)
