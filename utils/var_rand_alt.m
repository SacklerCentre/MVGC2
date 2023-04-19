function A = var_rand_alt(n,p,rho,lam,tol)

% Generate a random VAR coefficients sequence with given spectral radius
% and exponential decay factor.
%
% n       - observation variable dimension, or connectivity matrix/array
% p       - number of lags
% rho     - spectral radius
% lam     - exponential decay factor (empty for "natural" decay)
%
% A       - VAR coefficients sequence (3D array, last index is lag)

if nargin < 4,                 lam = [];        end
if nargin < 5 || isempty(tol), tol = sqrt(eps); end

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

if isempty(lam)
	A = specnorm(A,rho); % just enforce spectral norm
	return
end

% Enforce coefficients decay

for k = 1:p
	A(:,:,k) = (exp(-lam*k)/norm(A(:,:,k)))*A(:,:,k);
end

% Find appropriate weighting

wmin = 0;
wmax = 1;
% get upper bound on weight
while specnorm(wmax*A) < rho
	wmax = 2*wmax;
end
% binary chop
rhou = rho+tol/2;
rhol = rho-tol/2;
while true
	w = (wmin+wmax)/2;
	rhow = specnorm(w*A);
	if     rhow < rhol
		wmin = w;
	elseif rhow > rhou
		wmax = w;
	else
		break % within tolerance tol of rho
	end
end

% Strictly enforce spectral norm (more important than decay!)

A = specnorm(w*A,rho);
