n     = 10;     % total vars
p     = 6;    % model order

rho   = 0.99;   % spectral norm
w     = 1.5;   % decay weighting parameter (empty for no weighting)
g     = 0.5;   % residuals multi-information (0 for zero residuals correlation, empty for uniform random)

% Simulation parameters

N     = 10;     % number of trials
m     = 1000;   % number of observations (time series length)

%-------------------------------------------------------------------------------

% generate random VAR parameters satisfying H0

A0 = var_rand(n,p,rho,w); % VAR coefficients array
V0 = corr_rand(n,g);      % residuals covariance matrix

X = var_to_tsdata(A0,V0,m,N);

figure(1); clf
ptic('OLS\n');
tsdata_to_varmo(X,3*p,'OLS',[],[],[],2);
ptoc;

figure(2); clf
ptic('LWR\n');
tsdata_to_varmo(X,3*p,'LWR',[],[],[],2);
ptoc;
