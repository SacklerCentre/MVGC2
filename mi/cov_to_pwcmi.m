function I = cov_to_pwcmi(V)

% Pairwise-conditional mutual information between all pairs of variables,
% conditioned on remaining variables in the system.
%
% This quantity corresponds to the partial correlation between all pairs of
% variables in the system.
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% Distribution under the null hypothesis of zero MI for a particular pair is
% chi^2 with degrees of freedom d = 1, scaled by sample size = (number of
% trials) x (number of observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

LDVI = zeros(n,1);
for i = 1:n
    oi = 1:n; oi(i) = []; % omit i-th
    LDVI(i) = logdet(V(oi,oi));
end

LDVIJ = nan(n);
for i = 1:n
	for j = i+1:n
        oij = 1:n; oij([i j]) = []; % omit i-th and j-th
        LDVIJ(i,j) = logdet(V(oij,oij));
	end
end

I = symmetrise(LDVI+LDVI'-LDVIJ-logdet(V));
