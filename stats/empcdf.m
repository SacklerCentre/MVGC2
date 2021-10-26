function P = empcdf(stat,x)

% Quick-and-dirty empirical CDF. Calculates the empirical
% CDF of the statistic in each column of stat, evaluated at
% the points in the vector x.

assert(ismatrix(stat),'Statistics must be supplied as a matrix');
assert(isvector(x),'Evaluation points must be a vector');

n = size(stat,2);
m = length(x);
P = zeros(m,n);
x = x(:)';
for i = 1:n
	P(:,i) = mean(stat(:,i) <= x);
end
