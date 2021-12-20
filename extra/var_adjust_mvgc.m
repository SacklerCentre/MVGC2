function [AA,Fxy,iters] = var_adjust_mvgc(A,V,x,y,F,maxi,ftol);

% Find VAR coefficients AA based on supplied VAR parameters (A,V),
% such that Fxy, the conditional Granger causality from% y to x,
% is within a given tolerance ftol of the specified value F. The
% spectral radius of the model is preserved.
%
% The algorithm implements a binary chop. If the tolerance ftol
% is zero, the search continues until the binary chop can
% proceed no further, or the maximum number of iterations maxi
% is exceeded. Otherwise, the search continues until the tolerance
% is met, the binary chop cannot proceed, or the maximum number
% of iterations is exceeded. In either case, the returned result
% Fxy and number of iterations iters should be checked.

assert(F >= 0,'Target GC must be non-negative!');

% Reasonable defaults

if nargin < 6 || isempty(maxi), maxi = 1000;  end
if nargin < 7 || isempty(ftol), ftol = 1e-10; end

% Spectral radius

rho = specnorm(A);

% F == 0 is a special case:

if F == 0
	AA = A;
	AA(x,y,:) = 0;
	AA = specnorm(AA,rho);
	Fxy = 0;
	iters = 1;
	return
end

% GC(y -> x) is (generally!?) monotonic with scaling of the coefficients
% A(x,y,:); we try to find a multiplicative factor f which delivers the
% required GC. Firstly we find a ceiling for f.

Axy = A(x,y,:);
f = 0;
iters = 0;
while iters <= maxi
	f = f+1;
	AA = A;
	AA(x,y,:) = f*Axy;
	AA = specnorm(AA,rho);
	Fxy = var_to_mvgc(AA,V,x,y);
	if Fxy > F, break; end % have a ceiling for f
	iters = iters+1;
end

% We've boxed the multiplicative factor between f-1 and f; now binary chop
% as far as possible (i.e., to machine precision w.r.t. the factor f) or
% until Fxy is within the specified tolerance of F.

fmin = f-1;
fmax = f;
hftol = ftol/2;
while iters <= maxi
	fpre = f;                         % previous value of f
	f = (fmin+fmax)/2;                % binary chop
	if abs(fpre-f) <= eps, break; end % terminate if difference between old and new value of f is negligible
	AA = A;
	AA(x,y,:) = f*Axy;
	AA = specnorm(AA,rho);
	Fxy = var_to_mvgc(AA,V,x,y);
	if Fxy > F+hftol, fmax = f; elseif Fxy < F-hftol, fmin = f; else, break; end % terminate if tolerance met
	iters = iters+1;
end
