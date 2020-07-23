r   = 7;
n   = 5;
rho = 0.9;
g   = 0.5;
m   = 1000000;

pmax = 60;

q = 20;

%rng_seed(867817212);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A,C,K] = iss_rand(n,r,rho);

V = corr_rand(n,g);

ss_info(A,C,K,V);

X = ss_to_tsdata(A,C,K,V,m);
Y = X(:,m:-1:1);

Gx = tsdata_to_autocov(X,q);
Gy = tsdata_to_autocov(Y,q);
for k = 1:q+1
	Gyy(:,:,k) = Gy(:,:,k)';
end

Fx = autocov_to_pwcgc(Gx)
Fy = autocov_to_pwcgc(Gy)'

Fx-Fy
(Fx-Fy)./max(Fx+Fy)

return


[moaicx,mobicx,mohqcx,molrtx] = tsdata_to_varmo(X,pmax,'LWR',[],false,true,false,'');
[moaicy,mobicy,mohqcy,molrty] = tsdata_to_varmo(Y,pmax,'LWR',[],false,true,false,'');

p = round((molrtx+molrty)/2)

[Ax,Vx] = tsdata_to_var(X,p,'LWR');
[Ay,Vy] = tsdata_to_var(Y,p,'LWR');

Fx = var_to_pwcgc(Ax,Vx)
Fy = var_to_pwcgc(Ay,Vy)'

Fx-Fy
