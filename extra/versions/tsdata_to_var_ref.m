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
%     dm         de-mean data (default: true)
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
% *_Note_*: If the regressions are rank-deficient or ill-conditioned then A may
% be "bad" (i.e. will contain a |NaN| or |Inf|; see <isbad.html |isbad|>) and/or
% warnings will may be issued. The caller should test for both these
% possibilities. Possible causes are non-stationarity and/or colinearity in the
% data.
%
% The caller should also, at the very least, check the _spectral radius_ of the
% returned VAR coefficients (see <var_specrad |var_specrad|>) to ensure that the
% coefficients define a stable VAR [1]. (This is calculated, along with other
% relevant information, in the routine <var_to_autocov.html |var_to_autocov|>,
% which will typically be called subsequent to this function, and may be tested
% by a call to <var_info.html |var_info|>).
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
% <var_specrad.html |var_specrad|> |
% <var_to_autocov.html |var_to_autocov|> |
% <isbad.html |isbad|> |
% <var_info.html |var_info|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [A,V,E] = tsdata_to_var_ref(X,p,regmode,dm)

if nargin < 3 || isempty(regmode), regmode = 'OLS'; end
if nargin < 4 || isempty(dm),      dm      = true;  end

[n,m,N] = size(X);
assert(p < m,'too many lags');
p1 = p+1;

A = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')
V = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')
E = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')

if dm, X = demean(X); end % no constant term!

if  strcmpi(regmode,'OLS') % OLS (QR decomposition)

    M = N*(m-p);
    np = n*p;

    % stack lags

    X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
    XL = zeros(n,p,M);
    for k = 1:p
        XL(:,k,:) = reshape(X(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
    end
    XL = reshape(XL,np,M);         % stack lags

    A = X0/XL;                     % OLS (via QR decomposition)
    if isbad(A); return; end       % something went badly wrong

    if nargout > 1
        E = X0-A*XL;               % residuals
        V = (E*E')/(M-1);          % residuals covariance matrix (unbiased estimator)
        if nargout > 2             % align residuals per-trial with data (lose p lags)
			E = cat(2,nan(n,p,N),reshape(E,n,m-p,N));
		end
    end

    A = reshape(A,n,n,p);          % so A(:,:,k) is the k-lag coefficients matrix

elseif strcmpi(regmode,'LWR') % LWR (Morf)

    p1n = p1*n;

    I = eye(n);

    % store lags

    XX = zeros(n,p1,m+p,N);
    for k = 0:p
        XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
    end

    % initialise recursion - v2.0: EF/EB Initialisation corrected (many thanks to Gonzalo Camba-Mendez for the heads-up)

    EE = reshape(X,n,N*m);
    IC = inv(chol(EE*EE','lower')); % inverse covariance square root

    k  = 1;
    kn = k*n;
    M  = N*(m-k);
    kf = 1:kn;         % forward  indices
    kb = p1n-kn+1:p1n; % backward indices

    AF = zeros(n,p1n); % forward  AR coefficients
    AF(:,kf) = IC;

    AB = zeros(n,p1n); % backward AR coefficients (reversed compared with [2])
    AB(:,kb) = IC;

    % and loop

    while k <= p

        EF = AF(:,kf)*reshape(XX(:,1:k,k+1:m,:),kn,M); % forward  prediction errors
        EB = AB(:,kb)*reshape(XX(:,1:k,k:m-1,:),kn,M); % backward prediction errors

        CEF = chol(EF*EF');
        CEB = chol(EB*EB');

        R = CEF'\(EF*EB')/CEB; % normalised reflection coefficients

        [RF,cholp] = chol(I-R*R','lower');
        if cholp, return; end  % this should never happen!

        [RB,cholp] = chol(I-R'*R,'lower');
        if cholp, return; end  % this should never happen!

        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
        kf = 1:kn;
        kb = p1n-kn+1:p1n;

        AFPREV = AF(:,kf);
        ABPREV = AB(:,kb);

        AF(:,kf) = RF\(AFPREV-R*ABPREV);
        AB(:,kb) = RB\(ABPREV-R'*AFPREV);

    end

    % v2.0 - amended for an error that particularly affected p = 1. Note that
    % now, compared with OLS, we lose an extra lag at the front of E
    if nargout > 1
        E = AF(:,1:n)\AF(:,kf)*reshape(XX(:,1:k,k+1:m,:),kn,M); % residuals
		%E = reshape(X(:,p1:m,:),n,M) - A*reshape(XX(:,2:p1,p1:m,:),np,M); % alt - technically correct?
        V = (E*E')/(M-1); % residuals covariance matrix (unbiased estimator)
        if nargout > 2    % align residuals per-trial with data (lose p+1 lags)
			E = cat(2,nan(n,k,N),reshape(E,n,m-k,N));
		end
    end

    A = reshape(-AF(:,1:n)\AF(:,n+1:end),n,n,p); % so A(:,:,k) is the k-lag coefficients matrix

else
    error('bad regression mode ''%s''',regmode);
end
