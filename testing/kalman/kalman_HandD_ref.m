% Naive implementation of Kalman filter (Hannan & Deistler p. 90ff)

function [V,K,e,P,x] = kalman_HandD_ref(y,A,C,KK,VV)

[n,r] = ss_parms(A,C,KK,VV);

[n1,m] = size(y);
assert(n1 == n','Data doesn''t match SS parameters');

KL = KK*chol(VV,'lower');

Q = KL*KL';
S = KK*VV;
R = VV;

V = zeros(n,n,m);
K = zeros(r,n,m);
e = zeros(n,  m);
P = zeros(r,r,m);
x = zeros(r,  m); % x(t|t-1)

t = 1;
P(:,:,1) = dlyap(A,Q); % P_1 = A*P_1*A' + Q (assuming stationarity)
while true
	V(:,:,t) = C*P(:,:,t)*C' + R;
	K(:,:,t) = (A*P(:,:,t)*C'+S)/V(:,:,t);

	e(:,t) = y(:,t) - C*x(:,t);

	if t == m, break; end

	P(:,:,t+1) = A*P(:,:,t)*A' + Q - K(:,:,t)*V(:,:,t)*K(:,:,t)';
	x(:,t+1) =  A*x(:,t) + K(:,:,t)*e(:,t);

	t = t+1;
end
