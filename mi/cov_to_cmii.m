function I = cov_to_cmii(V,x)

% Multi-information is calculated for the multivariable specified by the
% vector x, conditioning on all other variables in the system.
%
% The covariance matrix V may be calculated parametrically from a VAR or SS
% model, or nonparametrically directly from the data.
%
% Distribution under the null hypothesis of zero MI is chi^2 with degrees of
% freedom d = nx*(nx-1)/2, scaled by sample size = (number of trials) x (number
% of observations per trial).

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

if nargin < 2 || isempty(x)
	x = 1:n; % all variables
end

x = x(:)'; % vectorise multivariable indices
assert(all(x >=1 & x <= n),'Some x indices out of range');
nx = length(x);

LDVOI = zeros(n,1);
ox = 1:n; ox(x) = [];
for i = x
	iox = [i ox];
	LDVOI(i) = logdet(V(iox,iox));
end
I = sum(LDVOI) - logdet(V) - (nx-1)*logdet(V(ox,ox));
