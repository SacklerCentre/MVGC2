function f = var_to_smvgc_ss(VARA,V,x,y,fres)


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
r = [x z];               % indices of reduced variables
xr = 1:length(x);        % index of x in reduced variables
yz = [y z];

[r,r1]  = size(A);   assert(r1 == r);
[n,r1]  = size(C);   assert(r1 == r);
[r1,n1] = size(K);   assert(n1 == n && r1 == r);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);
yz = [y z];

h = fres+1;
f = nan(1,h);

H   = ss2trfun(A,C,K,fres);
VL  = chol(V,'lower');
PVL = chol(parcov(V,yz,x),'lower');

% Note: in theory we shouldn't have to take the real part of the determinants of
% the (Hermitian, positive-definite) matrices in the calculation of the f(k),
% since they should be real. However, Matlab's det() will sometimes return a
% tiny imaginary part.

if isempty(z) % unconditional (note: does not require reduced model)

    for k = 1:h
        HVL  = H(:,:,k)*VL;
        SR   = HVL*HVL';
        HR   = H(x,y,k)*PVL;
        f(k) = logdet(SR(x,x)) - logdet(SR(x,x)-HR*HR');
    end

else % conditional

	% Solve the reduced DARE

    KVL = K*VL;
    CR  = C(r,:);
    [KR,VR,rep,~,P]  = ss2iss(A,CR,KVL*KVL',V(r,r),K*V(:,r)); % reduced model Kalman gain and innovations covariance
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

	% Calculate spectral GC

    BR    = ss2itrfun(A,CR,KR,fres);
    SR    = VR(xr,xr); % reduced model spectrum is flat!
    LDSR  = logdet(SR);
    for k = 1:h
        HR   = BR(xr,:,k)*H(r,yz,k)*PVL;
        f(k) = LDSR - logdet(SR-HR*HR');
    end

end
