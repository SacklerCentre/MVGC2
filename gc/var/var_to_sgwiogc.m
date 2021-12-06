function f = var_to_sgwiogc(A,V,group,inout,fres)

% In/out group spectralGCs
%
% Groups supplied as a cell vector of index vectors.
%
% For each group, GC is calculated either to or from the group and
% the rest of the system.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(group,n);

if strcmpi(inout,'in')
	gcin = true;
elseif strcmpi(inout,'out')
	gcin = false;
else
	error('in/out parameter must be ''in'' or ''out''');
end

h = fres+1;
f = nan(g,h);

H   = var2trfun(A,fres);
VL  = chol(V,'lower');

for a = 1:g
	if gcin
		x = group{a};
		y = 1:n; y(x) = [];
	else
		y = group{a};
		x = 1:n; x(y) = [];
	end
	PVL = chol(parcov(V,y,x),'lower');
    for k = 1:h
        HVL  = H(x,:,k)*VL;
        SR   = HVL*HVL';
        HR   = H(x,y,k)*PVL;
        f(a,k) = logdet(SR) - logdet(SR-HR*HR');
    end
end
