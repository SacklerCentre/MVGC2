function I = cov_to_gwcmii(V,groups)

% For each group, returns multi-information of group conditional on remaining
% variables in the system.
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% Distribution under the null hypothesis of zero MI for a group is chi^2 with
% degrees of freedom d = ng*(ng-1)/2, where ng is the size of the group, scaled
% by sample size = (number of trials) x (number of observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

[g,ng] = check_group(groups,n);

I = nan(g,1);
LDV  = logdet(V);
for a = 1:g
    goa = 1:n; goa(groups{a}) = [];
    I(a) = -LDV - (ng(a)-1)*logdet(V(goa,goa));
    for i = groups{a}
        igoa = [i goa];
        I(a) = I(a) + logdet(V(igoa,igoa));
    end
end
