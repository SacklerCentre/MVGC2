
% d = ???
% n = ???

[s,K] = LZc_maxcomp(d,n)

c = LZc(s,true,d)

if n < 0
	n1 = LZc_maxcompn(d,-n);
else
	n1 = n;
end
[length(s) n1 K]
