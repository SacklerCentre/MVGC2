function [Ak,C,Kk,Vk,rep] = ss_downsample_ver(A,C,K,V,k)

assert(isscalar(k) && isint(k) && k > 0,'Bad downsample factor (must be > 0)');

[n,r] = ss_parms(A,C,K,V);

U = K*chol(V,'lower');
KVK = U*U'; % ensure KVK = K*V*K' is really symm. pos-def.
Ak = eye(r);
Qk = KVK;
for m = 1:k-1
    Qk = A*Qk*A' + KVK;
    Ak = A*Ak;
end
Ak1KV = Ak*K*V; % Ak = A^(k-1)
Ak = A*Ak;      % Ak = A^k

[Kk,Vk,rep] = ss2iss(Ak,C,Qk,V,Ak1KV); % convert downsampled SS parms to ISS parms
