function F = var_to_gwcggc(A,V,groups)

% Groupwise-conditional GGCs
%
% Groups supplied as a cell vector of index vectors.
%
% Global GC is calculated for each group, conditioning on all other
% variables in the system specified by the VAR paremeters A,V.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(groups,n);

F = nan(g,1);

DV = diag(V);

for a = 1:g
	x = groups{a};
    z = 1:n; z(x) = [];
    nx = length(x);
	VRx = zeros(nx,1);
	for i = 1:nx
		y = x; y(i) = [];
		r = [x(i) z];
		[~,VR,rep] = vardare(A,V,y,r);
		if sserror(rep,i) % check DARE report, bail out on error
			sserr = true;
			break;
		end
		VRx(i) = VR(1,1);

	end
	F(a) = sum(log(VRx)) - sum(log(DV(x)));
end
