function I = cov_to_gwiomi(V,groups,m,N)

% For each group, returns mutual information between group and rest of system
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% Distribution under the null hypothesis of zero MI for a particular group
% is chi^2 with degrees of freedom d = ng*(n-ng), scaled by sample size
% = (number of trials) x (number of observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

g = check_group(groups,n);

LDV = logdet(V);

I = nan(g,1);
for a = 1:g
	ga = groups{a};
	goa = 1:n; goa(ga) = [];
	I(a) = logdet(V(ga,ga)) + logdet(V(goa,goa)) - LDV;
end
