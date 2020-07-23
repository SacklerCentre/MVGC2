% n = ???
% d = ???
% T = ???

usemex = 1;

for t = 1:T
	s = char(96+randi(d,1,n)); % high complexity
	c(:,t) = LZc_x(s);
end
cmean = mean(c,2);

cmax = LZc_cmaxx(n,d);

gp_qplot((1:n)',[cmean cmax]);
