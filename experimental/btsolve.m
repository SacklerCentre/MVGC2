% Block Levinson recursion for efficiently solving symmetric block-Toeplitz
% matrix equations.
%
% X = btsolve(G1,Y) solves the matrix equation T*X = Y, where T
% is a symmetric matrix with block-Toeplitz structure, and returns the
% solution matrix X. The matrix T is never stored in full (because it
% is large and mostly redundant); the input parameter G1 is the sequence
% of top block-row blocks.
%
% Adapted from 'block_levinson.m' by Keenan Pepper, 2007-12-23
%
% https://uk.mathworks.com/matlabcentral/fileexchange/30931-block-levinson-solver
%
% based on an algorithm described in:
%
% Hirotugu Akaike, "Block Toeplitz Matrix Inversion", SIAM J. Appl. Math. 24(2):234-241, 1973.

function X = btsolve(G1,Y)

[n,n1,q] = size(G1);     % n = block dimension, q = number of blocks
assert(n1 == n,'Blocks in Toeplitz input block sequence not square');
qn = q*n;

[qn1,m] = size(Y);
assert(qn1 == qn,'RHS input matrix does not match Toeplitz matrix');

G0 = G1(:,:,1);
L = G1(:,:);

G = reshape(permute(flipdim(G1,3),[2,1,3]),n,qn); % bottom block row of Toeplitz matrix

Id = eye(n);
Zn = zeros(n);
Zm = zeros(n,m);

X = zeros(qn,m);

% Initialise

F = inv(G0);             % "Forward" block vector
B = F;                   % "Backward" block vector
X(1:n,:) = F*Y(1:n,:);   % Solution vector

% Recurse

for k = 2:q
	kd = k*n;
	xidx = 1:kd;
	fidx = (k-1)*n+1:kd;
	bidx = (q-k)*n+1:qn;
	EF = G(:,bidx)*[F;Zn];
	EB = L(:,1:kd)*[Zn;B];
	I1 = inv(Id-EB*EF);
	I2 = inv(Id-EF*EB);
	AA = [[I1;-EF*I1],[-EB*I2;I2]];
	Fk = [[F;Zn],[Zn;B]]*AA(:,1:n);
	Bk = [[F;Zn],[Zn;B]]*AA(:,n+1:end);
	F  = Fk;
	B  = Bk;
	X(xidx,:) = X(xidx,:)+B*(Y(fidx,:)-G(:,bidx)*X(xidx,:));
end
