n = 1;
m = 10000;

X = randn(n,m);

Y0 = gendiff(X,[d r false]);
Y1 = gendiff(X,[d r true ]);

maxabs(Y0-Y1)

gp_qplot((1:m)',[Y0;Y1]',[],[],[],[1 Inf])
