% X = tsolve(G,Y) solves the matrix equation T*X = Y, where T
% is a symmetric Toeplitz matrix, and returns the solution matrix X.
% The input parameter G is the top row of T. Uses Levinson recursion.

function X = tsolve(G,Y)

assert(isrow(G),'Toeplitz matrix top row must be a row-vector');
q = length(G);

[q1,m] = size(Y);
assert(q1 == q,'Input matrix does not match Toeplitz row sequence');

G = flipdim(G,2);  % bottom row of Toeplitz matrix

X = zeros(q,m);

% Initialise

F = 1/G(end);      % forward  vector
B = F;             % backward vector
X(1,:) = F*Y(1,:); % solution matrix

% Recurse

for k = 2:q
	fidx = 1:k;
	bidx = q-k+1:q;
	F  = [F;0];
	B  = [0;B];
	E  = G(:,bidx)*F;
	M  = [F B]/(1-E*E);
	F  = M*[1;-E];
	B  = M*[-E;1];
	X(fidx,:) = X(fidx,:)+B*(Y(k,:)-G(:,bidx)*X(fidx,:));
end
