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

DV = diag(V);
LDV = log(DV);
I = LDV+LDV' - log(DV*DV'-V.*V);
I(1:n+1:n*n) = NaN; % make sure NaNs on diagonal.
