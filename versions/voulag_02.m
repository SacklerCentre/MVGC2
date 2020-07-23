function [X,t] = voulag(T,fs,NODE,CONX,V,usemex)

% Simulate lagged Ornstein-Uhlenbeck process X on a network, possibly with oscillatory nodes
%
%    dX(t) = -d X(t) dt + sum_(r = 1:c) a(r)*X(t-lag(r)) + dW(t)
%
% Z          n x m residuals sequence
%
% fs         simulation frequency (Hz)
%
% NODE       n x 2 matrix of node decay/oscillation frequencies (Hz)
%     fdec   decay frequencies (Hz)
%     fosc   oscillator frequencies (Hz)
%
% CONX       c x 4 matrix specifying connections, where c = number of connections: columns are
%     ito    "to"   node index, 1..n
%     ifr    "from" node index, 1..n
%     tlag   causal lag (seconds)
%     fcau   causal feedback frequency (strength) (Hz)
%
% usemex     use C implementation in 'oulag_mex.c (orders of magnitude faster!)
%
% The sample frequency fs should be "large" - say 10,000-100,000 Hz, and decay
% frequencies fdec << fs. For stability, we need fdec not "too small" and fcau
% not "too large" (absolute value).

if nargin < 6 || isempty(usemex), usemex = true; end

assert(isa(T,    'double') && isscalar(T)  && T  > 0, 'Simulation time must be a double-precision positive scalar');
assert(isa(fs,   'double') && isscalar(fs) && fs > 0, 'Sample rate must be a double-precision positive scalar');
assert(isa(NODE, 'double') && ismatrix(NODE),         'Decay/oscillator frequencies must be a double-precision matrix');
assert(isa(CONX, 'double') && ismatrix(CONX),         'Connectivity spec must be a double-precision matrix');
assert(isa(V,    'double') && ismatrix(V),            'Residuals covariance matrix must be a double-precision matrix');

[n,ncols] = size(NODE);
assert(ncols ==2,'Decay/oscillator frequencies mmatrix must have 2 columns');

dt = 1/fs;         % simulation time step
m = 1+round(T*fs); % number of observations

fdec = NODE(:,1);
fosc = NODE(:,2);

[c,ncols] = size(CONX);
assert(ncols == 4,'Connectivity matrix must have 4 columns');
ito  = CONX(:,1); % "to"   node index (1..n)
ifr  = CONX(:,2); % "from" node index (1..n)
tlag = CONX(:,3); % causal lag (SECONDS))
fcau = CONX(:,4); % causal feedback frequency (strength) (Hz)
uto = uint64(ito-1);
ufr = uint64(ifr-1);
assert(all([uto;ufr] < uint64(n)),'Some connectivity table nodes out of range');
ito = double(uto+1);
ifr = double(ufr+1);
assert(all(tlag > 0),'causal lags must be positive');

[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square and match nodes');
[L,cholp] = chol(V,'lower');
assert(cholp == 0,'Residuals covariance matrix must be positive-definite');

mlag = round(tlag*fs);

d = 1-2*pi*fdec*dt; % normalised decay frequencies (radians)
o = 2*pi*fosc*dt;   % normalised oscillator frequency (radians)
a = 2*pi*fcau*dt;   % normalised causal feedback frequencies

[d o]
a

ZX = L*randn(n,m)*sqrt(dt);
ZY = L*randn(n,m)*sqrt(dt);

if usemex
    X = voulag_mex(ZX,ZY,d,a,uto,ufr,uint64(mlag)); % implements precisely the code below
else
	X = ZX;
	Y = ZY;
	for t = 2:m
		X(:,t) = X(:,t) + d.*X(:,t-1) - o.*Y(:,t-1);
		Y(:,t) = Y(:,t) + o.*X(:,t-1) + d.*Y(:,t-1);
		for r = 1:c
			if mlag(r) < t
				X(ito(r),t) = X(ito(r),t) + a(r)*X(ifr(r),t-mlag(r));
			end
		end
	end
end

t = (0:m-1)*dt;
