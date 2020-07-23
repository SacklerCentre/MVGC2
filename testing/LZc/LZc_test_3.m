% d = ???
% n = ???

s0 = 'abcaaabacbabbbccacbccaaaaabaacabaabbabcacaacbaccbaababbacbbabbbbbcbcabcbbcccaacabcaccbacbbcbcccaccbccc';

s = s0;
s = 'abcdef';
%s = [s0 'a'];
%s = s(randperm(length(s)));

d = length(unique(s));
n = length(s);

d = 5; n = 100000; s = char(96+randi(d,1,n));
%d = 5; n = 100200; s = char(96+ones(1,n));

c  = LZc(s,true,d)

c1 = maxC1(n,d)
c2 = maxC(n,d)

floor((sqrt(8*n+1)-1)/2)

function c = maxC(n,d)

	g = 0; % g_k(d)     = d + 2*d^2 + ... + k*d^k
	f = 0; % f_{k-1}(d) = d + d^2 + ... + d^{k-1}
	dk = 1; % d^k
	k = 1;
	while true
		dk = dk*d;
		kdk = k*dk;
		if g+kdk > n % overshot: decrement k
			c = f + floor((n-gold)/k); %f_{k-1}(d) + D, D = floor((n-gold)/k)
			break
		end
		g = g+kdk;
		f = f+dk;
		k = k+1;
	end

end

function c = maxC1(n,d)

	k = 1;
	while true
		Dmax = realpow(d,k+1)-1;
		g = gkd(k,d);
		gotD = false;
		for D = 0:Dmax
			if n >= g + (k+1)*D && n < g+(k+1)*(D+1)
				gotD = true;
				break
			end
		end
		if gotD
			break
		end
		k = k+1;
	end
%[k D]
	c = fkd(k,d) + D;

end

function f = fkd(k,d)

	f = (d./(d-1))*(realpow(d,k)-1);

end

function g = gkd(k,d)

	d1 = d-1;
	g = (d./(d1*d1))*((k*d1-1)*realpow(d,k)+1);

end

%{
d = 5; n1 = 40; s1 = char(96+randi(d,1,n1));
d = 5; n2 = 55; s2 = char(96+randi(d,1,n2));
s = {s1;s2}
[c1,dict1] = LZc(s1,false)
%}
