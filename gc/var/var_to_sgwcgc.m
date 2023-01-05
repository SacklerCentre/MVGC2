function f = var_to_sgwcgc(A,V,groups,fres)

% Groupwise-conditional spectral GCs
%
% Groups supplied as a cell vector of index vectors.
%
% GC is calculated between each pair of groups (first index is target
% group, second is source group), conditioning on all other variables
% in the system specified by the VAR paremeters A,V.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(groups,n);

pn = p*n;
pn1 = pn-n;

h = fres+1;
f = nan(g,g,h);

H = var2trfun(A,fres);

% for efficiency (in time if not memory) we pre-compute partial covariances

PVL = cell(g,1);
for a = 1:g
    x = groups{a};
    w = 1:n; w(x) = []; % omit groups a
    PVL{a} = chol(parcov(V,w,x),'lower'); % pre-compute the partial covariances for efficiency
end

for b = 1:g
	y = groups{b};
    r = 1:n; r(y) = []; % omit groups b

	ny   = length(y);
	nr   = length(r);
	pny  = p*ny;
	pny1 = pny-ny;

	% Solve the shrunken DARE

	[KT,VR,rep] = vardare(A,V,y,r);
    if sserror(rep,b), continue; end % check DARE report, bail out on error

	% Calculate reduced SS parameters from shrunken DARE (note: VR is the same)

	AR = [reshape(A,n,pn); eye(pn1) zeros(pn1,n)];
	CR = reshape(A(r,:,:),nr,pn);

	%KR = kt2kr_mex(KT,r,y); % mex version of below: actually, no faster!

	KR = zeros(pn,nr);
	KR(r,:) = eye(nr);
	kn = 0;
	for ky = 0:ny:pny1
		KR(kn+y,:) = KT(ky+1:ky+ny,:);
		kn = kn+n;
	end
    BR = ss2itrfun(AR,CR,KR,fres);

    for a = 1:g
        if a == b, continue; end
        x = groups{a};
        xr = findin(x,r);   % indices of groups{a} in r
        w = 1:n; w(x) = []; % omit groups a

        SR  = VR(xr,xr);   % reduced model spectrum is flat!
        LDSR = logdet(SR);

        PVLa = PVL{a};
        for k = 1:h
            HR = BR(xr,:,k)*H(r,w,k)*PVLa;
            f(a,b,k) = LDSR - logdet(SR-HR*HR');
        end
    end
end
