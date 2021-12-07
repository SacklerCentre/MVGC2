function I = cov_to_mvmi(V,x,y)

% Mutual information between variables specified by indices x, y conditional
% on all other variables in the system.
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% Distribution under the null hypothesis of zero MI is chi^2 with degrees of
% freedom d = nx*ny, scaled by sample size = (number of trials) x (number of
% observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z = 1:n; z([x y]) = [];
xz = [x,z];
yz = [y,z];

I = logdet(V(xz,xz)) + logdet(V(yz,yz)) - logdet(V(z,z)) - logdet(V);
