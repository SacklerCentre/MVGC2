%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('n',   'var'), n   = 5;    end
if ~exist('r',   'var'), r   = 7;    end
if ~exist('rho', 'var'), rho = 0.9;  end
if ~exist('g',   'var'), g   = 0.4;  end
if ~exist('m',   'var'), m   = 1000; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A,C,KK] = iss_rand(n,r,rho);

VV = corr_rand(n,g);

y = ss_to_tsdata(A,C,KK,VV,m);

[V1,K1,e1,P1,x1] = kalman_HandD_ref(y,A,C,KK,VV);
[V,K,e,P,x] = kalman_HandD(y,A,C,KK,VV);
