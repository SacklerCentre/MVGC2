function [F,k,dF,tol] = ifgc_F_ver(A,C,K,V,x,y,stat,kmax,tol)

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

switch upper(stat)
	case 'LR', LRstat = true;
	case 'F',  LRstat = false;
	otherwise error('''stat'' must be ''LR'' or ''F''');
end

if nargin < 9 || isempty(tol), tol = sqrt(eps); end

F0 = NaN;

VL = chol(V,'lower');
KVL = K*VL;
[KR,VR,rep] = ss2iss(A,C(x,:),KVL*KVL',V(x,x),K*V(:,x));
if sserror(rep), return; end
CAk1 = C(x,:);

% Note: Dk = C(x,:)*A^{k-1}*K is the x component of the k-th MA coefficients matrix

F = nan(kmax,1);

if LRstat
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
else
	VRL   = chol(VR,'lower');
	KRVRL = KR*VRL;
	D     = VL(x,:); % D0
	DR    = VRL;
	F(1)  = trace(DR*DR')/trace(D*D')-1;
	for k = 2:kmax
		D    = [D  zeros((k-1)*nx,n ); CAk1*KVL   D( end-nx+1:end,:)];
		DR   = [DR zeros((k-1)*nx,nx); CAk1*KRVRL DR(end-nx+1:end,:)];
		Fk = sum(DR(:).^2)/sum(D(:).^2)-1;
		dF   = abs(F(k)-F(k-1));
		if dF < tol, break; end
		CAk1 = CAk1*A;
	end
%{
	VRL   = chol(VR,'lower');
	KRVRL = KR*VRL;
	Dk    = VL(x,:);
	t(1)  = trace(Dk*Dk');
	T(1)  = t(1);
	DRk   = VRL;
	tR(1) = trace(DRk*DRk');
	TR(1) = tR(1);
	F(1)  = TR(1)/T(1) - 1; % F(1) is 1-step GC
	for k = 2:kmax
		Dk    = CAk1*KVL;
		t(k)  = trace(Dk*Dk');
		T(k)  = sum((k:-1:1).*t(1:k));
		DRk   = CAk1*KRVRL;
		tR(k) = trace(DRk*DRk');
		TR(k) = sum((k:-1:1).*tR(1:k));
		F(k)  = TR(k)/T(k) - 1;
		dF    = abs(F(k)-F(k-1));
		if dF < tol, break; end
		CAk1 = CAk1*A;
	end
%}
end

F = F(1:k);
