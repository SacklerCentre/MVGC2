function [F,k,dF] = ss_to_ffgc_ver1(A,C,K,V,x,y,kmax,tol)

% Infinite future (total) GC of y -> x, up to kmax steps ahead (or until convergence)

% CAREFUL! D can grow to (kmax+1)^2 x nx x n

[n,r] = ss_parms(A,C,K,V);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

if nargin < 8 || isempty(tol), tol = -Inf; end % all 'kmax' values

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
nx = length(x);
xr = 1:nx;        % index of x in reduced quantities

VL = chol(V,'lower');
KVL = K*VL;
[KR,VR,rep] = ss2iss(A,C(x,:),KVL*KVL',V(x,x),K*V(:,x));
if sserror(rep), return; end

CAk1 = C(x,:);

% Note: Dk = C(x,:)*A^{k-1}*K is the x component of the k-th MA coefficients matrix

F = nan(kmax,1);

D = VL(x,:); % D0
LDVR = logdet(VR);
F(1) = LDVR-logdet(V(x,x)); % F(1) is 1-step GC
for k = 2:kmax
	D    = [D zeros((k-1)*nx,n); CAk1*KVL D(end-nx+1:end,:)];
	F(k) = k*LDVR-logdet(D*D');
	dF   = abs(F(k)-F(k-1));
	if dF < tol, break; end
	CAk1 = CAk1*A;
end

F = F(1:k);
