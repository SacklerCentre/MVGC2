function X = oulag(b,A,V,dt,lagt,simt,sett,usemex)

% X = oulag(b,A,V,dt,lagt,simt,sett, usemex = true)
%
% Simulate lagged Ornstein-Uhlenbeck process X
%
%    dX(t) = -b X(t) dt + sum_(r = 1:R) A(r)*X(t-lagt(r)) + dW(t)
%
% b          vector of decay parameters (must be > 0  for stability)
% A          square matrix of lagged coefficients
% V          residuals covariance matrix (must be square and positive-definite): dWTt) = V*sqrt(dt)
% dt         time increment (all times in seconds)
% lagt       lag times
% simt       simulation time
% sett       settle time
% usemex     use C implementation in 'oulag_mex.c (orders of magnitude faster!)

if nargin < 8 || isempty(usemex), usemex = true; end

[n,n1,R] = size(A);
assert(n1 == n,'coefficients matrices not square');
assert(isvector(lagt) && length(lagt) == R,'lags don''t match coefficient matrices');
if isscalar(b)
	b = b*ones(n,1);
else
	assert(isvector(b) && length(b) == n,'decay parameters must be a scalar or a vector matching the coefficients matrix');
end
assert(size(V,1) == n,'covariance matrix doesn''t match coefficients matrix');
assert(size(V,2) == n,'covariance matrix not square');

simm = round(simt/dt);
setm = round(sett/dt);
lagm = round(lagt/dt);
totm = simm+setm;
assert(all(diff(lagm) >= 1),'lags must be ascending (and not too close together!)');
blagm = lagm(R); % longest lag
assert(blagm < totm,'lag time(s) too large');

[C,p] = chol(V,'lower');
assert(p == 0,'covariance matrix not positive-definite');
Z = C*randn(n,totm)*sqrt(dt); % scaled Gaussian white noise

A = A*dt;
b = 1-b*dt;

if usemex
    X = oulag_mex(Z,b,A,lagm); % implements precisely the code below
else
    B = diag(b);
    X = Z;
    for t = blagm+1:totm
        Xt = X(:,t) + B*X(:,t-1);
        for r = 1:R
            Xt = Xt +A(:,:,r)*X(:,t-lagm(r));
        end
        X(:,t) = Xt;
    end
end

X = X(:,setm+1:totm);
