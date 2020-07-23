function f = ss_to_sgwcgc_alt(A,C,K,V,group,fres)

% inter-group (conditional) spectral GCs

[r,r1]  = size(A); assert(r1 == r);
[n,r1]  = size(C); assert(r1 == r);
[r1,n1] = size(K); assert(n1 == n && r1 == r);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

g = length(group);

h = fres+1;
f = nan(g,g,h);

H   = ss2trfun(A,C,K,fres);
KVL = K*chol(V,'lower');
KVK = KVL*KVL';

for b = 1:g
    r = 1:n; r(group{b}) = []; % omit group b

    CR = C(r,:);
    [KR,VR,rep] = ss2iss(A,CR,KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
    if sserror(rep,b), continue; end % check DARE report, bail out on error

    BR = ss2itrfun(A,CR,KR,fres);

    for a = 1:g
        if a == b, continue; end
        x = group{a};
        xr = findin(x,r);   % indices of group{a} in r
        w = 1:n; w(x) = []; % omit group a

        SR  = VR(xr,xr);   % reduced model spectrum is flat!
        LDSR = logdet(SR);

        PVL = chol(parcov(V,w,x),'lower'); % pre-compute the partial covariances for efficiency
        for k = 1:h
            HR = BR(xr,:,k)*H(r,w,k)*PVL;
            f(a,b,k) = LDSR - logdet(SR-HR*HR');
        end
    end
end
