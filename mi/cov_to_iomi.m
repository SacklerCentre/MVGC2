function I = cov_to_iomi(V)

% Mutual information between each variable and all remaining variables in
% the system.
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% For each variable, distribution under the null hypothesis of zero MI is chi^2
% with degrees of freedom d = n-1, scaled by sample size = (number of trials) x
% (number of observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

LDVOI = zeros(n,1);
for i = 1:n
	oi = 1:n; oi(i) = [];
	LDVOI(i) = logdet(V(oi,oi));
end
I = log(diag(V)) + LDVOI - logdet(V);
