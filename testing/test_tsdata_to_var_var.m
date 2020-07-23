n     = 5;     % total vars
p     = 1;    % model order

rho   = 0.999;   % spectral norm
w     = 1;   % decay weighting parameter (empty for no weighting)
g     = 0.5;   % residuals multi-information (0 for zero residuals correlation, empty for uniform random)

% Simulation parameters

N     = 3;     % number of trials
m     = 200;   % number of observations (time series length)

%rng_seed(214124);

S = 100;


%-------------------------------------------------------------------------------

% generate random VAR parameters satisfying H0


for s = 1:S

	fprintf('sample %d of %d\n',s,S);

	A0 = var_rand(n,p,rho,w); % VAR coefficients array
	V0 = corr_rand(n,g);      % residuals covariance matrix

	X = var_to_tsdata(A0,V0,m,N);
	D0 = det(V0);

	[A1,V1,E1] = tsdata_to_var(X,p,'LWR');
	fprintf('\tLWR NEW : rho = %6.4f  eA = %10.8f  eV = %10.8f\n',var_specrad(A1),maxabs(A1-A0),maxabs(V1-V0))
	D1(s) = det(V1)/D0;

	[A2,V2,E2] = tsdata_to_var_ref(X,p,'LWR');
	fprintf('\tLWR REF : rho = %6.4f  eA = %10.8f  eV = %10.8f\n',var_specrad(A2),maxabs(A2-A0),maxabs(V2-V0))
	D2(s) = det(V2)/D0;

	[A3,V3,E3] = tsdata_to_var(X,p,'OLS');
	fprintf('\tOLS     : rho = %6.4f  eA = %10.8f  eV = %10.8f\n',var_specrad(A3),maxabs(A3-A0),maxabs(V3-V0))
	D3(s) = det(V3)/D0;
end

mean(1-D1)
mean(1-D2)
mean(1-D3)

maxabs(A1-A2)
