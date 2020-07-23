function F = ss_to_pwcgc(A,C,K,V)

[n,~,L] = ss_parms(A,C,K,V);

F = nan(n);

KL  = K*L;
KVK = KL*KL';
LDV = log(diag(V));

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

    [~,VR,rep] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
    if sserror(rep,y), continue; end % check DARE report, bail out on error

    F(r,y) = log(diag(VR))-LDV(r);
end
