function [A,lam] = var_rand(n,p,rho,w,plotm)

% Generate a random VAR coefficients sequence with given spectral radius
% and VAR coefficients decay weighting factor.
%
% n       - observation variable dimension, or connectivity matrix/array
% p       - number of lags
% rho     - spectral radius
% w       - VAR coefficients decay weighting factor: empty (default) = don't weight
%
% A       - VAR coefficients sequence (3D array, last index is lag)
% lam     - VAR coefficients exponential decay factor

if nargin < 4, w     = []; end

if isscalar(n)
	A = randn(n,n,p);
else
	C = n; % connectivity matrix
	[n,n1,q] = size(C);
	assert(n1 == n,'Connectivity matrix must be square');
	if q == 1
		C = repmat(C,[1 1 p]);
	else
		assert(isempty(p),'If full connectivity array, model order must be empty');
		p = q;
	end
    A = C.*randn(n,n,p);
end

if isempty(w)
	A = specnorm(A,rho);
else
	A = specnorm(exp(-w*sqrt(p))*A,rho);
end

if nargout > 1
	lam = vardec(A);
end
