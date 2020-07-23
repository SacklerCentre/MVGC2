function [X,t] = voulag(T,fs,NODE,CONX,ONODE,OCONX,V,usemex)

% Simulate lagged Ornstein-Uhlenbeck process X on a network, possibly with oscillatory nodes
%
%    dX(t) = -d X(t) dt + sum_(k = 1:c) a(k)*X(t-lag(k)) + dW(t)
%
% Z          n x m residuals sequence
%
% fs         simulation frequency (Hz)
%
% NODE       n x 1 vector of node decay frequencies (Hz)
%     fdec   decay frequencies (Hz)
%
% CONX       c x 4 matrix specifying connections, where c = number of connections: columns are
%     ito    "to"   node index, 1..n
%     ifr    "from" node index, 1..n
%     fcau   feedback frequency (strength) (Hz)
%     tlag   causal lag (seconds)
%
% ONODE      on x 2 matrix of oscillator node decay/oscillation frequencies (Hz)
%     fodc   oscillator decay frequencies (Hz)
%     fosc   oscillator frequencies (Hz)
%
% OCONX      oc x 4 matrix specifying connections, where oc = number of oscillator connections: columns are
%     ioto   "to"   node index, 1..n
%     iofr   "from" oscillator node index, 1..on
%     foca   feedback frequency (strength) (Hz)
%
% usemex     use C implementation in 'oulag_mex.c (orders of magnitude faster!)
%
% The sample frequency fs should be "large" - say 10,000-100,000 Hz, and decay
% frequencies fdec << fs. For stability, we need fdec not "too small" and fcau
% not "too large" (absolute value).

if nargin < 6 || isempty(usemex), usemex = true; end

assert(isa(T, 'double') && isscalar(T)  && T  > 0, 'Simulation time must be a double-precision positive scalar');
assert(isa(fs,'double') && isscalar(fs) && fs > 0, 'Sample rate must be a double-precision positive scalar');

assert(isa(NODE,  'double') && isvector(NODE),'Node decay frequencies must be a double-precision vector');
assert(isa(CONX,  'double') && ismatrix(CONX),'Connectivity spec must be a double-precision matrix');
assert(isa(ONODE, 'double') && ismatrix(CONX),'Oscillator spec must be a double-precision matrix');
assert(isa(OCONX, 'double') && ismatrix(CONX),'Oscillator connectivity spec must be a double-precision matrix');

assert(isa(V,'double') && ismatrix(V),            'Residuals covariance matrix must be a double-precision matrix');

dt = 1/fs;         % simulation time step
m = 1+round(T*fs); % number of observations

n = length(NODE);
fdec = NODE(:);

[c,ncols] = size(CONX);
assert(ncols == 4,'Connectivity matrix must have 4 columns');
ito  = CONX(:,1); % "to"   node index (1..n)
ifr  = CONX(:,2); % "from" node index (1..n)
fcau = CONX(:,3); % causal feedback frequency (strength) (Hz)
tlag = CONX(:,4); % causal lag (SECONDS))
uto = uint64(ito-1);
ufr = uint64(ifr-1);
assert(all(uto < uint64(n)),'Some "to" nodes out of range');
assert(all(ufr < uint64(n)),'Some "from" nodes out of range');
ito = double(uto+1);
ifr = double(ufr+1);
assert(all(tlag > 0),'causal lags must be positive');
mlag = round(tlag*fs);

[on,ncols] = size(ONODE);
assert(ncols == 2,'Oscillator node matrix must have 2 columns');
fodc = ONODE(:,1);
fosc = ONODE(:,2);

[oc,ncols] = size(OCONX);
assert(ncols == 3,'Oscillator connectivity matrix must have 3 columns');
ioto = OCONX(:,1); % "to"   node index (1..n)
iofr = OCONX(:,2); % "from" node index (1..on)
foca = OCONX(:,3); % causal feedback frequency (strength) (Hz)
uoto = uint64(ioto-1);
uofr = uint64(iofr-1);
assert(all(uoto < uint64(n) ),'Some oscillator "to" nodes out of range');
assert(all(uofr < uint64(on)),'Some oscillator "from" nodes out of range');
ioto = double(uoto+1);
iofr = double(uofr+1);

[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square and match nodes');
[L,cholp] = chol(V,'lower');
assert(cholp == 0,'Residuals covariance matrix must be positive-definite');

d = 1-2*pi*fdec*dt; % normalised decay frequencies (radians)
a = 2*pi*fcau*dt;   % normalised causal feedback frequencies

od = 1-2*pi*fodc*dt; % normalised oscillator decay frequencies (radians)
oo = 2*pi*fosc*dt;   % normalised oscillator oscillation frequencies (radians)
oa = 2*pi*foca*dt;   % normalised oscillator causal feedback frequencies

ZX = L*randn(n,m)*sqrt(dt);
ZY = 1*randn(on,m)*sqrt(dt);
ZZ = 1*randn(on,m)*sqrt(dt);

if usemex
	error('sorry')
%    X = voulag_mex(ZX,ZY,d,a,uto,ufr,uint64(mlag)); % implements precisely the code below
else
	X = ZX;
	Y = ZY;
	Z = ZZ;
	for t = 2:m
		X(:,t) = X(:,t) + d.*X(:,t-1);
		Y(:,t) = Y(:,t) + od.*Y(:,t-1) - oo.*Z(:,t-1);
		Z(:,t) = Z(:,t) + oo.*Y(:,t-1) + od.*Z(:,t-1);
		for k = 1:oc
			X(ioto(k),t) = X(ioto(k),t) + oa(k)*Y(iofr(k),t-1);
		end
		for k = 1:c
			if mlag(k) < t
				X(ito(k),t) = X(ito(k),t) + a(k)*X(ifr(k),t-mlag(k));
			end
		end
	end
end

t = (0:m-1)*dt;
