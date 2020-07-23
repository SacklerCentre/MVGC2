function lags = decorrlags(rho,n,alpha)

% Heuristic autocorrelation decay length (based on Fisher
% z-transform of sample correlation statistic).
%
% Note: alpha smaller means critical value larger, so LESS lags required!
%
% rho   - autocorrelation decay factor (spectral radius)
% n     - number of sample observations
% alpha - correlation significance level

if nargin < 3 || isempty(alpha), alpha = 0.05; end

u = exp(-(2/sqrt(n-3))*norminv(alpha/2,0,1)); % for 2-tailed z-test
lags = ceil(log((u-1)/(u+1))/log(rho));
