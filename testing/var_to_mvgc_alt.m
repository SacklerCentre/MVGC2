function [F,rep] = var_to_mvgc_alt(VARA,V,x,y)

[n,n1,p] = size(VARA);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2]  = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
xz = [x z];
xr = 1:length(x);        % index of x in reduced quantities

ny = length(y);
nxz = length(xz);

pn1 = (p-1)*n;
pny = p*ny;
pny1 = pny-ny;

F = NaN;

A = [reshape(VARA(xz,,:),ny,pny); eye(pny1) zeros(pny1,ny)];
C = reshape(VARA(xz,y,:),nxz,pny);
Q = [V(y,y) zeros(ny,pny1); zeros(pny1,pny)];
S = [V(y,xz); zeros(pny1,nxz)];
R = V(xz,xz);
[K,VR,rep] = ss2iss(A,C,Q,R,S);
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
%CC = C;  CC(abs(CC)<eps) = 0
%K(abs(K)<eps) = 0
%HR = ss2trfun(A,C,K,100);
%S = HR(:,:,57)*VR*HR(:,:,57)'
F = logdet(VR(xr,xr)) - logdet(V(x,x));
