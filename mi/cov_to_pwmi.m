function I = cov_to_pwmi(V)

% Pairwise-unconditional mutual information between all pairs of variables.
%
% This quantity corresponds to the correlation between all pairs of
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

LDV = log(diag(V));
I = nan(n);
for i = 1:n
	for j = i+1:n
		ij = [i j];
		I(i,j) = LDV(i) + LDV(j) - logdet(V(ij,ij));
        I(j,i) = I(i,j);
	end
end
