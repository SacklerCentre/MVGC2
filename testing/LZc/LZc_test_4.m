% d = ???
% N = ???
% T = ???


c = zeros(N,1);
cx = zeros(N,1);
cc = zeros(T:1);
nn = (1:N)';
for n = 1:N
	fprintf('N = %d\n',n);

	for t = 1:T
		s = char(96+randi(d,1,n));
		cc(t) = LZc(s,true,d);
	end

	c(n) = mean(cc);

end

gp_qplot(nn,c);
