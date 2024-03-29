function f = ss_to_spwcgc(A,C,K,V,fres)

[n,r,L] = ss_parms(A,C,K,V);

h = fres+1;
f = nan(n,n,h);

H   = ss2trfun(A,C,K,fres);
KL  = K*L;
KVK = KL*KL';

% for efficiency (in time if not memory) we pre-compute partial covariances

nr = n-1;
PL = zeros(nr,nr,n);
for x = 1:n
    w = [1:x-1 x+1:n]; % omit x
    PL(:,:,x) = chol(parcov(V,w,x),'lower'); % pre-compute the partial covariances for efficiency
end

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

    CR = C(r,:);
    [KR,VR,rep] = mdare(A,CR,KVK,V(r,r),K*V(:,r)); % reduced model Kalman gain and innovations covariance
    if sserror(rep,y), continue; end % check DARE report, bail out on error

    BR = ss2itrfun(A,CR,KR,fres); % reduced inverse transfer function

    for xr = 1:n-1
        x  = r(xr);
        w = [1:x-1 x+1:n];  % omit x

        SR  = VR(xr,xr); % reduced model spectrum is flat!
        LSR = log(SR);
        PLx = PL(:,:,x);

        for k = 1:h
            HR = BR(xr,:,k)*H(r,w,k)*PLx;
            f(x,y,k) = LSR - log(SR-HR*HR');
        end
    end
end
