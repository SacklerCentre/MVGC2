function r = LZs(s,m)

% Lempl-Ziv-Welch relative complexity: relative, that is, to the mean
% complexity of a random string comprising exactly the same characters.
%
% INPUT
%
% s      input character string
% m      number of random samples
%
% OUTPUT
%
% r      LZW relative complexity

n  = LZW(s); % LZW complexity of string
l  = length(s);
nr = zeros(m,1);
for i = 1:m
	nr(i) = LZW(s(randperm(l))); % LZW complexity of randomly shuffled string
end

fprintf('relative deviation = %g\n',std(nr)/mean(nr));

r = n/mean(nr); % LZW relative complexity
