nx   = 3;
ny   = 5;
nz   = 2;
p    = 20;
rho  = 0.9;
w    = 1.0;
g    = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;
z = nx+ny+1:n;
pn1 = (p-1)*n;
r = [x z];
xr = x;        % index of x in reduced quantities

A0 = var_rand(n,p,rho,w);
V0 = corr_rand(n,g);

F0 = var_to_mvgc(A0,V0,x,y);

C = reshape(A0,n,p*n);
A = [C; eye(pn1) zeros(pn1,n)];
K = [eye(n); zeros(pn1,n)];

KVL = K*chol(V0,'lower');
[~,VR,rep] = ss2iss(A,C(r,:),KVL*KVL',V0(r,r),K*V0(:,r)); % reduced model innovations covariance
if sserror(rep), return; end % check DARE report, bail out on error
F1 = logdet(VR(xr,xr)) - logdet(V0(x,x));

A0R = A0; A0R(x,y,:) = zeros(nx,ny,p);
C2 = reshape(A0R,n,p*n);
A2 = [C2; eye(pn1) zeros(pn1,n)];

[~,VR2,rep] = ss2iss(A2,C2(r,:),KVL*KVL',V0(r,r),K*V0(:,r)); % reduced model innovations covariance
if sserror(rep), return; end % check DARE report, bail out on error
F2 = logdet(VR(xr,xr)) - logdet(V0(x,x));

[F0 F1 F2]
