function I = binary_mi(x,y)

% x, y should be matching binary vectors

n = length(x);
assert(isvector(x) && isvector(y) && length(y) == n, 'x and y must be binary vectors of the same length');

pxy = zeros(2,2); % probability table ("plugin" estimate)
for i = 0:1
	for j = 0:1
		pxy(i+1,j+1) = nnz(x == i & y == j)/n;
	end
end

px = sum(pxy'); % x marginal probabilities
py = sum(pxy);  % y marginal probabilities

I = entropy(px) + entropy(py) - entropy(pxy(:)); % mutual information

function H = entropy(p)

h = p.*log2(p);
h(p < eps) = 0; % p*log(p) --> 0 for p --> 0
H = -sum(h);
