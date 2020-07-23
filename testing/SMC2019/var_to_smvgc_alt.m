function f = var_to_smvgc_alt(VARA,V,x,y,fres)

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

nx = length(x);
ny = length(y);
nz = length(z);
nr = nx+nz;
xr = 1:nx;

r = [x z];
w = [y z];

pn = p*n;
pn1 = pn-n;
pny = p*ny;
pny1 = pny-ny;
pnr  = p*nr;

h = fres+1;
f = nan(1,h);

H   = var2trfun(VARA,fres);
VL  = chol(V,'lower');
PVL = chol(parcov(V,w,x),'lower');

if isempty(z) % unconditional (note: does not require reduced model)

    for k = 1:h
        HVL  = H(:,:,k)*VL;
        SR   = HVL*HVL';
        HR   = H(x,y,k)*PVL;
        f(k) = logdet(SR(x,x)) - logdet(SR(x,x)-HR*HR');
    end

else % conditional

	% Solve the shrunken DARE

	AT = [reshape(VARA(y,y,:),ny,pny); eye(pny1) zeros(pny1,ny)];
	CT = reshape(VARA(r,y,:),nr,pny);
	QT = [V(y,y) zeros(ny,pny1); zeros(pny1,pny)];
	ST = [V(y,r); zeros(pny1,nr)];
	RT = V(r,r);
	[KT,VR,rep,~,PT] = ss2iss(AT,CT,QT,RT,ST);
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

	% Calculate reduced DARE from shrunken DARE (note: VR is the same)

	AR = [reshape(VARA,n,pn); eye(pn1) zeros(pn1,n)];
	CR = reshape(VARA(r,:,:),nr,pn);

	MT = zeros(n,ny); MT(y,:) = eye(ny);
	NT = zeros(n,nr); NT(r,:) = eye(nr);
	KR = kron(eye(p),MT)*KT;
	KR(1:n,:) = KR(1:n,:) + NT;

	% Calculate spectral GC

    BR    = ss2itrfun(AR,CR,KR,fres);
    SR    = VR(xr,xr); % reduced model spectrum is flat!
    LDSR  = logdet(SR);
    for k = 1:h
        HR   = BR(xr,:,k)*H(r,w,k)*PVL;
        f(k) = LDSR - logdet(SR-HR*HR');
    end

end

%-------------------------------------------------------

%{
nr
ny
p
size(KT)
size(KR)

	KR1 = kt2kr_mex(KT,p);
size(KR1)
KR
KR1
%}
