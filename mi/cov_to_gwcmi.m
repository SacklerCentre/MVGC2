function I = cov_to_gwcmi(V,groups)

% For each pair of groups, returns mutual information between them,
% conditional on remaining variables in the system
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% Distribution under the null hypothesis of zero MI for a particular pair of
% groups is chi^2 with degrees of freedom d = nga*ngb, scaled by sample size
% = (number of trials) x (number of observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

g = check_group(groups,n);

LDV = logdet(V);

LDVG = zeros(g,1);
for a = 1:g
    goa = 1:n; goa(groups{a}) = [];
    LDVG(a) = logdet(V(goa,goa));
end

I = nan(g);
for a = 1:g
	for b = a+1:g
        goab = 1:n; goab([groups{a} groups{b}]) = [];
        I(a,b) = LDVG(a) + LDVG(b) - logdet(V(goab,goab)) - LDV;
	end
end
I = symmetrise(I);
