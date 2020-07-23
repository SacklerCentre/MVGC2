function [S,H] = ss_to_cpsd(A,C,K,V,fres)

[n,r,L] = ss_parms(A,C,K,V);

h  = fres+1;
In = eye(n);
Ir = eye(r);
S  = zeros(n,n,h);
w  = exp(1i*pi*((0:fres)/fres));
if nargout > 1
	H = zeros(n,n,h);
	for k = 1:h % over [0,pi]
		H(:,:,k) = In + C*((w(k)*Ir-A)\K);
		HLk      = H(:,:,k)*L;
		S(:,:,k) = HLk*HLk';
	end
else
	for k = 1:h % over [0,pi]
		HLk      = (In + C*((w(k)*Ir-A)\K))*L;
		S(:,:,k) = HLk*HLk';
	end
end
