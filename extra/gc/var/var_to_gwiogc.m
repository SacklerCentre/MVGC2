function F = var_to_gwiogc(A,V,groups,inout)

% In/out group GCs
%
% Groups supplied as a cell vector of index vectors.
%
% For each group, GC is calculated either to or from the group and
% the rest of the system.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(groups,n);

if strcmpi(inout,'in')
	gcin = true;
elseif strcmpi(inout,'out')
	gcin = false;
else
	error('in/out parameter must be ''in'' or ''out''');
end

F = nan(g,1);

for a = 1:g
	if gcin
		x = groups{a};
		y = 1:n; y(x) = [];
	else
		y = groups{a};
		x = 1:n; x(y) = [];
	end
	[~,VR,rep] = vardare(A,V,y,x);
	if sserror(rep), return; end % check DARE report, bail out on error
	F(a) = logdet(VR) - logdet(V(x,x));
end
