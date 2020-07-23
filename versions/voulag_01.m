function [X,t,Z] = voulag(NODE,CON,V,dt,simt,sett,usemex)

% Simulate lagged Ornstein-Uhlenbeck process X on a network, possibly with oscillatory nodes
%
%    dX(t) = -b X(t) dt + sum_(r = 1:c) a(r)*X(t-lag(r)) + dW(t)
%
% NODE       n x 2 matrix of decay/oscillation parameters, where n = number of nodes:
%            the first column contains the decay parameter (MILLISECONDS), the second
%            column contains the oscillation frequency (Hz); if the oscillation frequency
%            is non-zero, then a "virtual node" is used to generate a sinusoidal component
%
% CON        c x 4 matrix specifying connections, where c = number of connections: columns are
%     to     "to"   node
%     from   "from" node
%     lagt   causal lag (MILLISECONDS)
%     cfbs   causal feedback strength (1/SECONDS)
%
% V          n x n residuals covariance matrix (positive-definite): dWt) = V*sqrt(dt)
% dt         simulation time-step size (MILLISECONDS)
% simt       simulation time (SECONDS)
% sett       transient settle time (SECONDS)
% usemex     use C implementation in 'oulag_mex.c (orders of magnitude faster!)
%
% The time slice dt should be "small" - say 0.01 - 0.001, and NODE >> dt.
% For stability, we need NODE and cfbs not too large. Ballpark magnitudes:
%
% dt   ~ 0.01
% NODE ~ 40
% lagt ~ whatever
% cfbs ~ 20

if nargin < 7 || isempty(usemex), usemex = true; end

assert(ismatrix(NODE),'Oscillation/decay spec must be a matrix');
assert(ismatrix(CON), 'Connectivity spec must be a matrix');
assert(ismatrix(V),   'Residuals covariance matrix must be a matrix');

[n,ncols] = size(NODE);
assert(ncols == 2,'Oscillation/decay matrix must have 2 columns');
dect = NODE(:,1)/1000;        % node decay times (SECONDS)
oscf = NODE(:,2);             % node oscillator frequencies (Hz)
oidx = find(abs(oscf) > eps); % indices of oscillator nodes
nosc = length(oidx);          % number of oscillator nodes
ntot = n+nosc;                % total number of nodes, including virtual oscillator nodes

[c,ncols] = size(CON);
assert(ncols == 4,'Connectivity matrix must have 4 columns');
to   = CON(:,1);      % "to"   node
from = CON(:,2);      % "from" node
lagt = CON(:,3)/1000; % causal lag (SECONDS))
cfbs = CON(:,4);      % causal feedback strength (Hz)
assert(all([to;from] > 0) && all([to;from] <= n), 'Some connectivity table nodes out of range');

dt = dt/1000; % simulation step time (SECONDS)

[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match number of nodes');
V = [V zeros(n,nosc); zeros(nosc,n) eye(nosc)]; % duplicate for virtual oscillator nodes
[L,cholp] = chol(V,'lower');
assert(cholp == 0,'Residuals covariance matrix must be positive-definite');

simm = round(simt/dt)+1;
setm = round(sett/dt);
totm = simm+setm;

lagm = round(lagt/dt);

Z = L*randn(ntot,totm)*sqrt(dt); % scaled Gaussian white noise; the last nosc rows are for virtual oscillator nodes

d = 1-dt./[dect;dect(oidx)];
o = 2*pi*oscf*dt
a = dt*cfbs;

if usemex
    X = voulag_mex(Z,d,o,a,uint64(oidx-1),uint64(to-1),uint64(from-1),uint64(lagm)); % implements precisely the code below
else
oidx
	X = Z;
	for t = 2:totm
		X(:,t) = X(:,t) + d.*X(:,t-1); % decay all nodes, including virtual
		for k = 1:nosc                 % generate sinusoidal components for oscillator nodes
			i = oidx(k);
			j = n+k; % virtual node indices
			X(i,t) = X(i,t) - o(k)*X(j,t-1);
			X(j,t) = X(j,t) + o(k)*X(i,t-1);
		end
		for r = 1:c
			if lagm(r) < t
				X(to(r),t) = X(to(r),t) + a(r)*X(from(r),t-lagm(r));
			end
		end
	end
end

% truncate virtual oscillator nodes, and transients
X = X(1:n,setm+1:totm);
t = (0:simm-1)*dt;
