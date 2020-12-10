% X = btsolve(G,Y) solves the matrix equation T*X = Y, where T
% is a symmetric matrix with block-Toeplitz structure, and returns the
% solution matrix X. The matrix T is never stored in full (because it
% is large and mostly redundant); the input parameter G is the sequence
% (3d array) of top block-row blocks.
%
% Adapted from 'block_levinson.m' by Keenan Pepper, 2007-12-23
%
% https://uk.mathworks.com/matlabcentral/fileexchange/30931-block-levinson-solver
%
% based on a block-Levinson recursion algorithm described in:
%
% Hirotugu Akaike, "Block Toeplitz Matrix Inversion", SIAM J. Appl. Math. 24(2):234-241, 1973.

function X = btsolve(G,Y)

[n,n1,q] = size(G);      % n = block dimension, q = number of blocks
assert(n1 == n,'Blocks in Toeplitz block sequence not square');
qn = q*n;

[qn1,m] = size(Y);
assert(qn1 == qn,'Input matrix does not match Toeplitz matrix');

GF = G(:,:);                                      % top    block-row of Toeplitz matrix
GB = reshape(permute(flipdim(G,3),[2,1,3]),n,qn); % bottom block-row of Toeplitz matrix

In = eye(n);
Zn = zeros(n);
X  = zeros(qn,m);

% Initialise

F = inv(G(:,:,1));       % forward  block-vector
B = F;                   % backward block-vector
X(1:n,:) = F*Y(1:n,:);   % solution matrix

% Recurse

for k = 2:q
	kn = k*n;
	yidx = kn-n+1:kn;
	fidx = 1:kn;
	bidx = qn-kn+1:qn;
	F  = [F;Zn];
	B  = [Zn;B];
	FB = [F B];
	EF = GB(:,bidx)*F;
	EB = GF(:,fidx)*B;
	F  = FB*([In;-EF]/(In-EB*EF));
	B  = FB*([-EB;In]/(In-EF*EB));
	X(fidx,:) = X(fidx,:)+B*(Y(yidx,:)-GB(:,bidx)*X(fidx,:));
end
