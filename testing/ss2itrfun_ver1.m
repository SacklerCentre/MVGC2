function J = ss2itrfun_ver1(A,C,K,fres)

disp('ss2itrfun_v1')

p = 10000;

[n,r] = size(C);
h = fres+1;
B = A-K*C;
VARA = zeros(n,n,p);
Bk1 = eye(r);
for k = 1:p % over [0,pi]
	VARA(:,:,k) = C*Bk1*K;
	BK1 = Bk1*B;
end

J = fft(cat(3,eye(n),VARA),2*fres,3); % over [0,2*pi)
J = J(:,:,1:(fres+1));              % over [0,pi] only
