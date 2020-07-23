function r = gencorr(V,ascc)

% Generalised correlation is calculated as r = sqrt(1-exp(-g)), where g = -2log|R|/n
% with R the correlation matrix and n the number of dimensions; i.e.
%
% r = sqrt(1-|R|^(2/n))
%
% For n = 2 and small r, this is close to the Pearson correlation coefficient.
%
% For g, set ascc flag to false.

if nargin < 2 || isempty(ascc), ascc = true; end

r = 2*(sum(log(diag(V)))-logdet(V))/size(V,1);

if ascc
	r = sqrt(1-exp(-r));
end
