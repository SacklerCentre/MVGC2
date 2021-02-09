% Solve (OLS) the p-lag cross-regression problem
%
% Y_t = A_1 X_{t-1} + A_2 X_{t-2} + ... + A_p X_{t-p} + E_t
%
% for a jointly-stochastic process (X,Y) from the joint
% autocovariance sequence.
%
% INPUT:
%
%    G is the autocovariance sequence to p lags (see
%    'var_to_autocov', etc.), while x,y are index vectors
%    for the X,Y sub-processes of the joint process.
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
p  = size(G,3)-1;

Gxx = G(x,x,1:end-1);
Gyx = G(y,x,2:end);
Gyx = Gyx(:,:)';

if nx == 1 % Gxx is Toeplitz, as opposed to block-Toeplitz - use faster algorithm

	A = tsolve(Gxx(:)',Gyx)';

else

	A = btsolve(Gxx,Gyx)';

end

if nargout > 1
	V = G(y,y,1)-A*Gyx;
end

A = reshape(A,ny,nx,p);
