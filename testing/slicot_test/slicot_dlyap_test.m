% n = 8;

%-------------------------------------------------------------


A = randn(n);

C = symmetrise(randn(n));

tic
X = dlyap(A,-C);
et = toc; fprintf('MATLAB dlyap: %g secs\n',et);


tic
X1 = slicot_dlyap(A,C);
et = toc; fprintf('SLICOT dlyap: %g secs\n',et);

fprintf('|| A*X*A''-X -C ||   = %e\n',maxabs(A*X*A'-X -C));
fprintf('|| A*X1*A''-X1 -C || = %e\n',maxabs(A*X1*A'-X1 -C));
fprintf('|| X1-X ||          = %e\n',maxabs(X1-X));
