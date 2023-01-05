%% tsdata_to_var
%
% Fit VAR model to multi-trial, multivariate time series data
%
% <matlab:open('tsdata_to_var.m') code>
%
%% Syntax
%
%     [A,V,E] = tsdata_to_var(X,p,regmode)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     p          model order (number of lags)
%     regmode    regression mode: 'OLS' (default) or 'LWR'
%
% _output_
%
%     A          VAR coefficients matrix
%     V        residuals covariance matrix
%     E          residuals time series
%
%% Description
%
% Returns VAR coefficients |A| and (optionally) residuals covariance matrix
% |V| and serially uncorrelated residuals |E| for the |p|-lag autoregression
%
% <<eq_var.png>>
%
% (where  [[ii_Sigma.png]] = |V|) of a stationary multivariate process
% |X|. |X| may contain single- or multi-trial multivariate time series
% data. The regression mode is set by the |regmode| parameter, which may be
% |'LWR'| (default) or |'OLS'|. The former uses Morf's version of the LWR
% algorithm [1,2] while the latter calculates the OLS solution to the
% regression via QR decomposition.
%
% *_Note_*: If the regressions are rank-deficient or ill-conditioned then |A| may
% contain NaNs or Infs; this may be checked with the <isbad.html |isbad|>
% function. Possible causes are non-stationarity and/or colinearity in the data.
%
% The OLS mode is actually valid for unstable VAR processes, in particular unit
% root processes. The caller may check the _spectral radius_ of the returned VAR
% coefficients (see <specnorm.html |specnorm|>) to test whether the coefficients
% define a stable VAR [1]. (This is calculated, along with other relevant
% information, by a call to <var_info.html |var_info|>).
%
% For the LWR mode, the spectral radius is guaranteed to be < 1; this is because
% the algorithm _assumes_ stability. This implies that the result is invalid if
% the generative process is unstable (e.g.,a unit root process).
%
% Many thanks to Gonzalo Camba-Mendez for pointing out an initialisation bug in
% the original version (MVGC1) of the LWR code.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] M. Morf, A. Viera, D. T. L. Lee and T. Kailath, "Recursive Multichannel
% Maximum Entropy Spectral Estimation", _IEEE Trans. Geosci. Elec._, 16(2), 1978.
%
%% See also
%
% <specnorm.html |specnorm|> |
% <var_to_autocov.html |var_to_autocov|> |
% <var_info.html |var_info|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [A,V,E] = tsdata_to_var(X,p,regmode)

[n,m,N] = size(X);
assert(p < m,'too many lags or bad model order (p = %d, m = %d)',p,m);

M = N*(m-p); % effective number of observations

if p == 0 % white noise!
	A = zeros(n,n,0);
	if nargout > 1
		X0 = reshape(demean(X),n,M);
		V = (X0*X0')/(M-1);
		if nargout > 2
			E = X; % the residuals are just the time series itself!
		end
	end
	return
end

X = demean(X); % no constant term

if  strcmpi(regmode,'OLS') % OLS (QR decomposition)

	obs = p+1:m;  % the useable observations (i.e., lose the first p)

	X0 = reshape(X(:,obs,:),n,M); % concatenate trials for unlagged observations
	XL = zeros(n,p,M);
	for k = 1:p % for each lags
		XL(:,k,:) = reshape(X(:,obs-k,:),n,M); % concatenate trials for k-lagged observations
	end
	XL = reshape(XL,p*n,M); % stack lagged observations

	A = X0/XL; % OLS (via QR decomposition)

	if nargout > 1
		E = X0-A*XL; % residuals
	end

	A = reshape(A,n,n,p); % so A(:,:,k) is the k-lag coefficients matrix

elseif strcmpi(regmode,'LWR') % LWR (Morf et al.)

	I = eye(n);

	p1  = p+1;
	p1n = p1*n;

	% store lags

	XX = zeros(n,p1,m+p,N);
	for k = 0:p
		XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
	end

	% initialise recursion

	EE = reshape(X,n,N*m);

	[C,cholp] = chol(EE*EE','lower');
	assert(cholp == 0,'Covariance matrix not positive-definite (is there colinearity in your data?)');
	IC = inv(C); % inverse covariance square root

	k  = 1;
	kn = k*n;
	M  = N*(m-k);
	kk = 1:k;
	kf = 1:kn;         % forward  indices
	kb = p1n-kn+1:p1n; % backward indices

	AF = zeros(n,p1n); AF(:,kf) = IC; % forward  AR coefficients
	AB = zeros(n,p1n); AB(:,kb) = IC; % backward AR coefficients (reversed compared with [2])

	% LWR recursion

	while k <= p

		EF = AF(:,kf)*reshape(XX(:,kk,k+1:m,:),kn,M); % forward  prediction errors
		EB = AB(:,kb)*reshape(XX(:,kk,k:m-1,:),kn,M); % backward prediction errors

		R = (chol(EF*EF','lower')\EF)*(chol(EB*EB','lower')\EB)'; % normalised reflection coefficients

		k  = k+1;
		kn = k*n;
		M  = N*(m-k);
		kk = 1:k;
		kf = 1:kn;         % forward  indices
		kb = p1n-kn+1:p1n; % backward indices

		AFPREV = AF(:,kf);
		ABPREV = AB(:,kb);

		AF(:,kf) = chol(I-R*R','lower')\(AFPREV-R*ABPREV);
		AB(:,kb) = chol(I-R'*R,'lower')\(ABPREV-R'*AFPREV);

	end

	A0 = AF(:,1:n);
	A = reshape(-A0\AF(:,n+1:p1n),n,n,p); % so A(:,:,k) is the k-lag coefficients matrix

	if nargout > 1
		E = A0\EF;        % residuals
	end

else
	error('unknown regression mode ''%s''',regmode);
end

if nargout > 1
	V = (E*E')/(M-1); % residuals covariance matrix (unbiased estimator)
	if nargout > 2
%		E = cat(2,zero(n,p,N),reshape(E,n,m-p,N)); % pad first p lags of residuals with zeros to align with X
		E = reshape(E,n,m-p,N); % temp fix for whiteness, rsquared and consistency routines - better to pad!
	end
end
