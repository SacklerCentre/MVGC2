% n = 8;
% N = 100;

%-------------------------------------------------------------


A = randn(n,n,N);
C = randn(n,n,N);
for i = 1:N
	C(:,:,i) = symmetrise(C(:,:,i));
end

X  = zeros(n,n,N);
X1 = zeros(n,n,N);

tic
for i = 1:N
	X(:,:,i) = dlyap(A(:,:,i),-C(:,:,i));
end
et = toc; fprintf('MATLAB dlyap: %g secs\n',et);


tic
for i = 1:N
	X1(:,:,i) = slicot_dlyap(A(:,:,i),C(:,:,i));
end
et = toc; fprintf('SLICOT dlyap: %g secs\n',et);

W  = zeros(n,n,N);
W1 = zeros(n,n,N);
for i = 1:N
	W(:,:,i) = A(:,:,i)*X(:,:,i)*A(:,:,i)'-X(:,:,i) - C(:,:,i);
end
for i = 1:N
	W1(:,:,i) = A(:,:,i)*X1(:,:,i)*A(:,:,i)'-X1(:,:,i) - C(:,:,i);
end

fprintf('|| A*X*A''-X -C ||   = %e\n',maxabs(W));
fprintf('|| A*X1*A''-X1 -C || = %e\n',maxabs(W1));
fprintf('|| X1-X ||          = %e\n',maxabs(X1-X));
