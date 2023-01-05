function V = tsdata_to_cov(X,unbiased)

if nargin < 2 || isempty(unbiased), unbiased = true; end

X = X(:,:); % concatenate trials
X  = bsxfun(@minus,X,mean(X,2));
if unbiased
	V = (X*X')/(size(X,2)-1);
else
	V = (X*X')/size(X,2);
end
