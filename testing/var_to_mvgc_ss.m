function F = var_to_mvgc_ss(VARA,V,x,y)

[n,n1,p] = size(VARA);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2]  = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

[A,C,K] = var_to_ss(VARA,V);

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);        % index of x in reduced quantities

F = NaN;

KVL = K*chol(V,'lower');

[K,VR,rep,L,P] = ss2iss(A,C(r,:),KVL*KVL',V(r,r),K*V(:,r)); % reduced model innovations covariance
if rep < 0 % show-stopper!
    fprintf(2,'ERROR in reduced model calculation: ');
    switch rep
        case -1, fprintf(2,'DARE eigenvalues on/near unit circle\n');
        case -2, fprintf(2,'couldn''t find stablising DARE solution\n');
    end
    return
end
if rep > sqrt(eps)
    fprintf(2,'WARNING in reduced model calculation: DARE is inaccurate (relative residual = %e)\n',rep);
end

F = logdet(VR(xr,xr)) - logdet(V(x,x));

%-------------------------------------------------------

%{
ny = length(y);
nr = length(r);

M = zeros(ny,n); M(:,y) = eye(ny);
KT = kron(eye(p),M)*K

NT = zeros(n,nr); NT(r,:) = eye(nr);
K1 = kron(eye(p),M')*KT;
K1(1:n,:) = K1(1:n,:) + NT;

zdisp(K)
zdisp(K1)

maxabs(K-K1)
%}

%{
zdisp(P)
M = zeros(ny,n); M(:,y) = eye(ny);
MT = zeros(n,ny); MT(y,:) = eye(ny);
isequal(M',MT)
U = kron(eye(p),M);
UT = kron(eye(p),MT);
isequal(U',UT)
M
U
PT = U*P*U';
PT
zdisp(U'*PT*U)
N = zeros(nr,n); N(:,r) = eye(nr);
W = kron(eye(p),N);
N
W
J = [W;U];
J
zdisp(W*P*W')
zdisp(J*P*J')
zdisp(A)
zdisp(U*A*U')
zdisp(W*A*W')
zdisp(PT)
zdisp(W*A*U')
zdisp(W*A*U'*PT*(W*A*U')')
Q = KVL*KVL';
zdisp(W*Q*W')

y
v = zeros(n,1); v(y,1) = ones(ny,1)
MM = diag(v)
isequal(MM,M'*M)

UU = kron(eye(p),MM);

KT = U*K;



K1 = U'*KT;
K1(1:n,1:nr) = K1(1:n,1:nr) + N';
zdisp(K)
zdisp(K1)
maxabs(K-K1)
%}

%{
KKT = reshape(KT',[nr ny p]);
zdisp(KT)
for k = 1:p, disp(zdisp(KKT(:,:,k)')); end
%}

%{
KK = zeros(p*n,nr);
KK(1:n,:) = N';
KK() = KT;
zdisp(KK)
%}
