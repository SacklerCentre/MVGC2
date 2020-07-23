% Naive implementation of Kalman filter (Hannan & Deistler p. 90ff)

function [V,K,e,P,x] = kalman_HandD(y,A,C,KK,VV)

[n,r] = ss_parms(A,C,KK,VV);

[n1,m] = size(y);
assert(n1 == n','Data doesn''t match SS parameters');

[LL,p] = chol(VV,'lower');
assert(p == 0,'t = 0 : not posdef');

KKLL = KK*LL;

Q = KKLL*KKLL';
S = KK*VV;
R = VV;

L = zeros(n,n,m);
K = zeros(r,n,m);
e = zeros(n,  m);
U = zeros(r,r,m);
x = zeros(r,  m); % x(t|t-1)

t = 1;
[U(:,:,1),p] = chol(dlyap(A,Q),'lower'); % P_1 = A*P_1*A' + Q (assuming stationarity)
assert(p == 0,'t = %d : not posdef',t);
while true
	P(:,:,t) = U(:,:,t)*U(:,:,t)';

	CU = C*U(:,:,t);
	AU = A*U(:,:,t);

	[L(:,:,t),p] = chol(CU*CU' + R,'lower');
	assert(p == 0,'t = %d : not posdef',t);
	V(:,:,t) = L(:,:,t)*L(:,:,t)';

	KL = (AU*CU'+S)/L(:,:,t)';
	K(:,:,t) = KL/L(:,:,t);

	e(:,t) = y(:,t) - C*x(:,t);

	if t == m, break; end

	[U(:,:,t+1),p] = chol(AU*AU' + Q - KL*KL','lower');
	assert(p == 0,'t = %d : not posdef',t);

	x(:,t+1) =  A*x(:,t) + K(:,:,t)*e(:,t);

	t = t+1;
end
