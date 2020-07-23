function [F,k,dF,tol] = ifgc_trace(A,C,K,V,x,y,total,kmax,tol)

% CAREFUL! D can grow to (kmax+1)^2 x nx x n

[n,r] = ss_parms(A,C,K,V);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);        % index of x in reduced quantities

nx = length(x);

if nargin < 9 || isempty(tol), tol = sqrt(eps); end

F0 = NaN;

VL = chol(V,'lower');
KVL = K*VL;
[KR,VR,rep] = ss2iss(A,C(x,:),KVL*KVL',V(x,x),K*V(:,x));
if sserror(rep), return; end

CAk1 = C(x,:);

% Note: Dk = C(x,:)*A^{k-1}*K is the x component of the k-th MA coefficients matrix

VRL   = chol(VR,'lower');
KRVRL = KR*VRL;

T  = zeros(kmax,1);
TR = zeros(kmax,1);
Dk  = VL(x,:);
DRk = VRL;
T(1)  = trace(Dk*Dk');
TR(1) = trace(DRk*DRk');
for k = 2:kmax
	Dk    = CAk1*KVL;
	DRk   = CAk1*KRVRL;
	T(k)  = T(k-1)  + trace(Dk*Dk');
	TR(k) = TR(k-1) + trace(DRk*DRk');
	CAk1  = CAk1*A;
end

F = nan(kmax,1);

F(1) = (TR(1)-T(1))/T(1);


if total
	Sk  = T(1);
	SRk = TR(1);
	for k = 2:kmax
		Sk    = Sk  + T(k);
		SRk   = SRk + TR(k);
		Fk    = (SRk-Sk)/Sk;
		F(k)  = Fk;
		dF    = abs(F(k)-F(k-1))/(abs(F(k))+abs(F(k-1)));
		if dF < tol, break; end
	end
else
	for k = 2:kmax
		Fk    = (TR(k)-T(k))/T(k);
		F(k)  = F(k-1) + Fk;
		dF    = abs(Fk);
		if dF < tol, break; end
		CAk1 = CAk1*A;
	end
end

F = F(1:k);
