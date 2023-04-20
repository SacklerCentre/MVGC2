function [A,lrerr,evals] = var_rand_alt(n,p,rho,lam,tol)

% Generate a random VAR coefficients sequence with given spectral radius
% and exponential decay factor.
%
% n       - observation variable dimension, or connectivity matrix/array
% p       - number of lags
% rho     - spectral radius
% lam     - VAR coefficients exponential decay factor (empty for "natural" decay)
% tol     - binary chop tolerance
%
% A       - VAR coefficients sequence (3D array, last index is lag)
% lrerr   - relative error of actual decay factor lam
% evals   - binary chop spectral norm/decay evaluations

if nargin < 4,                  lam  = [];        end
if nargin < 5 || isempty(tol),  tol  = sqrt(eps); end

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

% Enforce spectral radius

A = specnorm(A,rho); % just enforce spectral norm

if isempty(lam) % just enforce spectral norm
	lrerr = [];
	evals = [];
	return
end

% Find appropriate weighting

% get upper bound on weight
w = 1;
Aw = specnorm(w*A,rho);
lamw = vardec(w*A);
evals = 1;
while lamw < lam
	w = 2*w;
	Aw = specnorm(w*A,rho);
	lamw = vardec(Aw);
	evals = evals+1;
end

% binary chop
wmin = 0;
wmax = w;
while true
	w = (wmin+wmax)/2;
	Aw = specnorm(w*A,rho);
	evals = evals+1;
	lamw = vardec(Aw);
	if     lamw < lam
		wmin = w;
	elseif lamw > lam
		wmax = w;
	else
		break % within tolerance tol of lam
	end
	if wmax-wmin < tol
		break % run out of wiggle-room!
	end
end
A = Aw;
lrerr = abs(lamw-lam)/abs(lam);
