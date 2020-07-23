function SAP = cpsd2apsd(SCP)

% Extract auto-power SD from cross-power SD

[n,n1,h] = size(SCP);
assert(n1 == n,'bad dimensions');

N = n*n;
SAP = SCP((1:n+1:N)+(0:h-1)'*N); % all the diagonals!

%{
SAP = zeros(h,n);
for k = 1:h
	SAP(k,:) = diag(SCP(:,:,k));
end
%}
