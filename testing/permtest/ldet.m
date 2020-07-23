function [D,flag] = ldet(V)

[L,flag] = chol(V);
if flag == 0 % symmetric, positive-definite
	D = prod(diag(L))^2;
else
	D = det(V);
end
