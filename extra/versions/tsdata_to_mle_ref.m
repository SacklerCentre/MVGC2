function [LL,Lk,Lm] = tsdata_to_mle_ref(X,q,regmode,verb)

[n,m,N] = size(X);

assert(isscalar(q) && isint(q) && q > 0 && q < m,'maximum model order must be a positive integer less than the number of observations');

X = demean(X); % no constant term

% store lags

q1 = q+1;
XX = zeros(n,q1,m+q,N);
for k = 0:q
    XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
end

% Note: k = 1 is order 0 term!
LL = nan(q1,1); % log-likelihood
Lk = nan(q1,1); % number free parameters
Lm = nan(q1,1); % effective sample size

if  strcmpi(regmode,'OLS') % OLS

    % order zero likelihood

    M  = N*m;
    X0 = reshape(X,n,M);        % concatenate trials for unlagged observations
    DSIG = det((X0*X0')/(M-1)); % covariance matrix determinant
    assert(DSIG > 0,'covariance matrix not positive-definite');

    LL(1) = -(M/2)*log(DSIG);
    Lk(1) =  0;
    Lm(1) =  M;

    % loop through model orders

    for k = 1:q

        if verb, fprintf('model order = %d',k); end

        k1 = k+1;
        M  = N*(m-k);

		X0 = reshape(XX(:,1,   k1:m,:),n,  M);
		XL = reshape(XX(:,2:k1,k1:m,:),n*k,M);

        wstate = warning('off','all'); lastwarn('');
        A = X0/XL;                     % OLS using QR decomposition
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(A)                     % something went badly wrong
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            continue % show-stopper
        end

        E    = X0-A*XL;                % residuals
        DSIG = det((E*E')/(M-1));      % residuals covariance matrix determinant

        if DSIG <= 0
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: residuals covariance not positive definite\n');
            continue % show-stopper
        end

        LL(k1) = -(M/2)*log(DSIG);
        Lk(k1) =  k*n*n;
        Lm(k1) =  M;

        if verb, fprintf(1,'\n'); end
    end

elseif strcmpi(regmode,'FLS') % Fast OLS (single QR decomposition) % v2.0 - new fast algorithm, based on QR decomposition

    % perform the QR decomposition for all lags in one shot

    M  = N*(m-q);
    X0 = reshape(XX(:,1,   q1:m,:),n,  M);
    XL = reshape(XX(:,2:q1,q1:m,:),n*q,M);
    [Q,R] = qr(XL',0);
    XQ = X0*Q;

    % order zero likelihood

    DSIG = det((X0*X0')/(M-1)); % covariance matrix determinant
    assert(DSIG > 0,'covariance matrix not positive-definite');

    LL(1) = -(M/2)*log(DSIG);
    Lk(1) =  0;
    Lm(1) =  M;

    % loop through model orders

    for k = 1:q

        if verb, fprintf('model order = %d',k); end

        k1 = k+1;
        nk = n*k;
        r = min(nk,M);
        wstate = warning('off','all'); lastwarn('');
        A = XQ(:,1:r)/R(1:r,1:nk)';     % OLS for regression against first k lags only
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(A)                     % something went badly wrong
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            continue % show-stopper
        end

        E    = X0-A*XL(1:nk,:);        % residuals
        DSIG = det((E*E')/(M-1));      % residuals covariance matrix determinant
        if DSIG <= 0
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: residuals covariance not positive definite\n');
            continue % show-stopper
        end

        LL(k1) = -(M/2)*log(DSIG);
        Lk(k1) =  k*n*n;
        Lm(k1) =  M;

        if verb, fprintf(1,'\n'); end
    end

elseif strcmpi(regmode,'LWR') % LWR (Levinson, Wiggins & Robinson algorithm - Morf variant)

    % order zero likelihood

    M  = N*m;
    DSIG = det((X(:,:)*X(:,:)')/(M-1)); % covariance matrix determinant
    assert(DSIG > 0,'covariance matrix not positive-definite');

    LL(1) = -(M/2)*log(DSIG);
    Lk(1) =  0;
    Lm(1) =  M;

    q1n = q1*n;

    I = eye(n);

    % initialise recursion

    AF = zeros(n,q1n); % forward  AR coefficients
    AB = zeros(n,q1n); % backward AR coefficients (reversed compared with Morf's treatment)

    k  = 1;            % model order is k-1
    kn = k*n;
    M  = N*(m-k);
    kf = 1:kn;         % forward  indices
    kb = q1n-kn+1:q1n; % backward indices

    EF = reshape(XX(:,1:k,k+1:m,:),kn,M);
    EB = reshape(XX(:,1:k,k:m-1,:),kn,M);

    [CEF,cholp] = chol(EF*EF');
    assert(cholp == 0,'initialisation failed (''forward'' covariance matrix not positive definite)'); % v2.0 - it's a show-stopper!

    [CEB,cholp] = chol(EB*EB');
    assert(cholp == 0,'initialisation failed (''backward'' covariance matrix not positive definite)'); % v2.0 - it's a show-stopper!

    wstate = warning('off','all'); lastwarn('');
    AF(:,kf) = CEF'\I;
    AB(:,kb) = CEB'\I;
    wmsg = lastwarn; warning(wstate);
    assert(isempty(wmsg),'initialisation failed (%s)\n',wmsg); % v2.0 - it's a show-stopper!

    % and loop

    while k <= q

        if verb, fprintf('model order = %d',k); end

        EF = AF(:,kf)*reshape(XX(:,1:k,k+1:m,:),kn,M); % forward  prediction errors
        EB = AB(:,kb)*reshape(XX(:,1:k,k:m-1,:),kn,M); % backward prediction errors

        [CEF,cholp] = chol(EF*EF');
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        [CEB,cholp] = chol(EB*EB');
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        R = CEF'\(EF*EB')/CEB;       % normalised reflection coefficients

        [CRF,cholp] = chol(I-R*R');
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        [CRB,cholp] = chol(I-R'*R);
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
        kf = 1:kn;
        kb = q1n-kn+1:q1n;

        AFPREV = AF(:,kf);
        ABPREV = AB(:,kb);

        wstate = warning('off','all'); lastwarn('');
        AF(:,kf) = CRF'\(AFPREV-R*ABPREV);
        AB(:,kb) = CRB'\(ABPREV-R'*AFPREV);
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg)
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed (%s)\n',wmsg);
            break % it's a show-stopper!
        end

        wstate = warning('off','all'); lastwarn('');
        E = AF(:,1:n)\AF(:,kf)*reshape(XX(:,1:k,k+1:m,:),kn,M);
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg)
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed (%s)\n',wmsg);
            break % it's a show-stopper!
        end

        DSIG = det((E*E')/(M-1));

        i = k-1;
        if DSIG <= 0
            if ~verb, fprintf('model order = %d',i); end
            fprintf(2,'  WARNING: residuals covariance matrix not positive definite\n');
            break % show-stopper
        end

        LL(i+1) = -(M/2)*log(DSIG);
        Lk(i+1) =  i*n*n;
        Lm(i+1) =  M;

        if verb, fprintf(1,'\n'); end
    end

else
    error('bad regression mode ''%s''',regmode);
end
