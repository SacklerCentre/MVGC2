function [I,pval] = cov_to_pwcmi(V,m)

% Pairwise-conditional MIs (partial correlation!)
%
% NOTE: if multi-trial, m = nobs x ntrials

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

LDV = logdet(V);

LDVI = zeros(n,1);
for i = 1:n
    oi = 1:n; oi(i) = []; % omit i-th
    LDVI(i) = logdet(V(oi,oi));
end

I = nan(n);

for i = 1:n
	for j = i+1:n
        oij = 1:n; oij([i j]) = []; % omit i-th and j-th
        I(i,j) = LDVI(i) + LDVI(j) - logdet(V(oij,oij)) - LDV;
        I(j,i) = I(i,j);
	end
end

if nargout > 1;
    pval  = 1-chi2cdf(m*I,1); % dof = 1
end
