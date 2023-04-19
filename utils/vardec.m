function [A,a] = vardec(A,lam)

% VAR coefficients sequence day
%
% [lam,a] = vardec(A) returns the fitted exponential decay factor
%           for the given VAR coefficients sequence A.
%
% [A,a]   = vardec(A,lam) returns the coefficient sequence A weighted
%           so that the decay factor = lam.
%
% Optionally, the norms at each lag are returned in a.
%
% Note that the second calling form will NOT in general preserve the
% spectral radius of the coefficients sequence.

p = size(A,3);
a = zeros(p,1);
for k = 1:p
	a(k) = norm(A(:,:,k));
end

if nargin < 2 || isempty(lam) % no decay factor supplied; return fitted decay factor

	parms = [(1:p)' ones(p,1)]\log(a); % fit to exponential (linear regression = OLS)
	A = -parms(1);                     % this is the fitted decay factor lam

else                          % weight the AR coefficients so decay = lam

	for k = 1:p
		A(:,:,k) = (exp(-lam*k)/a(k))*A(:,:,k);
	end

end
