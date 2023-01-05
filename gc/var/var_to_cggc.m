function F = var_to_cggc(A,V,x)

% Conditional GGC
%
% GGC is calculated for the multivariable specified by the vector x,
% conditioning on all other variables in the system specified by the
% VAR paremeters A,V.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

if nargin < 3 || isempty(x)
	x = 1:n; % all variables
end

x = x(:)'; % vectorise multivariable indices
assert(all(x >=1 & x <= n),'Some x indices out of range');

DV = diag(V);
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
F = sum(log(VRx)) - sum(log(DV(x)));
