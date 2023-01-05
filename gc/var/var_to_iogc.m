function F = var_to_iogc(A,V,inout)

% In/out GCs per variable
%
% For each variable, GC is calculated either to or from the variable
% and the rest of the system.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

if strcmpi(inout,'in')
	gcin = true;
elseif strcmpi(inout,'out')
	gcin = false;
else
	error('in/out parameter must be ''in'' or ''out''');
end

F = nan(n,1);
if gcin
	for x = 1:n
		y = 1:n; y(x) = [];
		[~,VR,rep] = vardare(A,V,y,x);
		if sserror(rep), return; end % check DARE report, bail out on error
		F(i) = log(VR) - log(V(x,x));
	end
else
	for y = 1:n
		x = 1:n; x(y) = [];
		[~,VR,rep] = vardare(A,V,y,x);
		if sserror(rep), return; end % check DARE report, bail out on error
		F(i) = logdet(VR) - logdet(V(x,x));
	end
end
