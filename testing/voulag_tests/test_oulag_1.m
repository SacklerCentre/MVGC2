n      = 1000;           % number of variables
nlags  = 3;           % number of lags

dtm    = 0.01;        % OU time-slice increment (MILLISECONDS)

lagtm  = [10 17 22];  % causal feedback lag times (MILLISECONDS)
mdectm = 4;           % mean decay time (MILLISECONDS)
sdectm = 2;           % std. dev. decay time (MILLISECONDS)

simt   = 30;          % simulation time (SECONDS)
sett   = simt;        % transient decay time (SECONDS)

usemex = true;       % use C (mex) version -- much faster

%-------------------------------------------------------------

assert(length(lagtm) == nlags);

% convert millisecond times to seconds

dt    = dtm/1000;
lagt  = lagtm/1000;
mdect = mdectm/1000;
sdect = sdectm/1000;

% generate Gamma-distributed decay times

ga = (mdect^2)/(sdect^2); % Gamma shape parameter
gb = (sdect^2)/mdect;     % Gamma scale parameter
dect = gamrnd(ga,gb,n,1); % Gamma-distributed decay-times


X = oulag(1./dect,A,V,dt,lagt,simt,sett,usemex);
