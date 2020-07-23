% n = ???
% d = ???
% T = ???

nn = (1:n)';
dd = (1:D)';

c = zeros(n,D,R);
for r = 1:R
	fprintf('sample %3d of %3d\n',r,R);
	x = randn(n,1);
	for d = 1:D
		s = LZc_quantise(x,d);
		c(:,d,r) = LZc_x(s);
	end
end
cmean = mean(c,3);

cmin = repmat(LZc_cmin(nn),1,D);

cmax = zeros(n,D);
for d = 1:D
	cmax(:,d) = LZc_cmaxx(n,d);
end

cnorm = (cmean-cmin)./(cmax-cmin);

gp_qplot((1:n)',cmean,num2str(dd,'d = %d'),'set title "mean LZc"\nset key bottom right\nset logs xy','png',1.5,16);

gp_qplot((1:n)',cmax,num2str(dd,'d = %d'),'set title "max LZc"\nset key bottom right\nset logs xy','png',1.5,16);

gp_qplot((1:n)',cnorm,num2str(dd,'d = %d'),'set title "normalised LZc"\nset key bottom right\nset logs x','png',1.5,16);
