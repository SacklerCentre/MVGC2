% Solve (OLS) the p-lag regression problem
%
% Y_t = A_1 X_{t-1} + A_2 X_{t-2} + ... + A_p X_{t-p} + E_t
%
% for a jointly-stochastic process (X,Y) from the joint
% autocovariance sequence.
%
% INPUT:
%
%    G is the autocovariance sequence to p lags (see
%    'var_to_autocov', etc.), while x,y are the index vectors
%    to the X,Y sub-processes of the joint process.
%
% OUTPUT:
%
%    A is the sequence of regression coefficients, V is the
%    residuals covariance matrix.
%
% The OLS solution is calculated in the 'btsolve' routine, which
% uses an algorithm due to H. Akaike for solving block-Toeplitz
% matrix equations.

function [A,V] = cross_regress(G,x,y)

nx = length(x);
ny = length(y);
p1 = size(G,3);
p = p1-1;

G0yy = G(y,y,1);
G1yx = G(y,x,2:p1);
G1xx = G(x,x,1:p);

Lxy = G1yx(:,:)';

A = btsolve(G1xx,Lxy)';

if nargout > 1
	V = G0yy-A*Lxy;
end

A = reshape(A,ny,nx,p);
