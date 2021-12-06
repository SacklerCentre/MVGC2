function f = var_to_scggc(A,V,x)

% Conditional spectral GGC
%
% GGC is calculated for the multivariable specified by the vector x,
% conditioning on all other variables in the system specified by the
% VAR paremeters A,V.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise multivariable indices
assert(all(x >=1 & x <= n),'Some x indices out of range');

TODO
