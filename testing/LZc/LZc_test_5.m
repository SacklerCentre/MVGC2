% N = ???


for n = 1:N
	fprintf('N = %d\n',n);

	s = char(97+zeros(1,n));

	c(n) = LZc(s,false);
	c1(n) = floor(sqrt(8*n+1)/2);
	c2(n) = floor((sqrt(8*n+1)-1)/2);

end

[c' c1' c2']
