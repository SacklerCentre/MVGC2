function B = var2vma_1(A,q)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');

B = zeros(n,n,q);

B(:,:,1) = A(:,:,1);
for k = 2:q
	C = zeros(n);
	for m = max(1,k-p):(k-1)
		C = C + A(:,:,k-m)*B(:,:,m);
	end
	B(:,:,k) = A(:,:,k) + C;
end
