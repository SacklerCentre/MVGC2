function [AA,FF,iters] = var_adjust_mvgc(A,V,x,y,F,maxi);

% Find VAR coefficients based on supplied VAR parameters (A,V)
% such that GC(y -> x) is equal to specified value F, based on
% a binary chop. Terminate when binary chop cannot proceed any
% further, or maximum iterations maxi is exceeded. Result should
% be tested; that is, check that FF lies as close as required to
% target FT.

assert(F >= 0,'Target GC must be non-negative!');

if nargin < 6 || isempty(maxi), maxi = 1000;  end

% Initialise

rho = specnorm(A);
iters = 0;

% F = 0 is a special case:

if F < eps
	AA = A;
	AA(x,y,:) = 0;
	AA = specnorm(AA,rho);
	FF = 0;
	return
end

% GC(y -> x) is (generally!?) monotonic with scaling of the coefficients
% A(x,y,:); we try to find a multiplicative factor f which delivers the
% required GC. Firstly we find a ceiling for f.

Axy = A(x,y,:);
f = 0;
FF = 0;
while FF < F
	iters = iters+1;
	if iters > maxi, return; end % timed out
	f = f+1;
	AA = A;
	AA(x,y,:) = f*Axy;
	AA = specnorm(AA,rho);
	FF = var_to_mvgc(AA,V,x,y);
end
fmin = f-1;
fmax = f;

% We've boxed f between fmin and fmax; now binary chop as far as possible
% (i.e., to machine precision w.r.t. f) or as far as necessary (i.e., FF is
% within machine precision of F).

htol = eps/2;
while true
	fpre = f;                         % previous value of f
	f = (fmin+fmax)/2;                % binary chop
	if abs(fpre-f) <= eps, break; end % no point going any further!
	iters = iters+1;
	if iters > maxi, return; end % timed out
	AA = A;
	AA(x,y,:) = f*Axy;
	AA = specnorm(AA,rho);
	FF = var_to_mvgc(AA,V,x,y);
	if FF > F+htol, fmax = f; elseif FF < F-htol, fmin = f; else, break; end % new box
end
