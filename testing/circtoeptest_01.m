
% a0 = ???
% b  = ???
% n = ???

bb = b^2;

a = zeros(n,1);
a(1) = a0-bb/a0;
for k = 2:n
	a(k) = a(k-1)-bb/a(k-1);
end

gp_qplot((0:n)',[a0;a]);
