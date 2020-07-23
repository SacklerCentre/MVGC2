function [F0,F,k,dF,tol] = ifgc(A,C,K,V,x,y,maxm,tol)

% CAREFUL! D can grow to (maxm+1)^2 x nx x n

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

if nargin < 8 || isempty(tol), tol = sqrt(eps); end

F0 = NaN;

VL = chol(V,'lower');
KVL = K*VL;
[~,VR,rep] = ss2iss(A,C(x,:),KVL*KVL',V(x,x),K*V(:,x));
if sserror(rep), return; end
LDVR = logdet(VR);

D = VL(x,:); % D0
CAk1 = C(x,:);
F0 = LDVR-logdet(V(x,x)); % 1-step GC
F = nan(maxm,1);
for k = 1:maxm
	D = [D zeros(k*nx,n); CAk1*KVL D(end-nx+1:end,:)]; % Bk = C*A^{k-1}*K is the k-th MA coefficients matrix
	F(k) = (k+1)*LDVR-logdet(D*D');
	if k > 1
		dF = F(k)-F(k-1);
		if dF < tol
			break;
		end
	end
	CAk1 = CAk1*A;
end
