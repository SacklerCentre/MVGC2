function [s,K] = LZc_maxcomp(n,d)

A = char(96+(1:d));

g = 0;
D = 1;
s = A;
if n < 0
	K = -n;
	for k = 1:K-1
		D = d*D;
		h = g;
		g = h+k*D;
		i = h;
		for m = 1:D
			sm = s(i+1:i+k);
			for j = 1:d
				s = [s sm A(j)];
			end
			i = i+k;
		end
	end
else
	N = d;
	k = 1;
	while true
		D = d*D;
		h = g;
		g = h+k*D;
		i = h;
		for m = 1:D
			sm = s(i+1:i+k);
			for j = 1:d
				s = [s sm A(j)];
				N = N+k+1;
				if N >= n
					s = s(1:n);
					K = k+1;
					return
				end
			end
			i = i+k;
		end
		k = k+1;
	end
end
