r   = 7;
nx  = 5;
ny  = 6;
rho = 0.95;
g   = 1;

m   = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

[A,C,K] = iss_rand(n,r,rho);

[V,C0] = corr_rand(n,g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VL = chol(V,'lower');
KVL = K*VL;
D0 = VL(x,:); % D0

% Initialise

DDk = D0;
dk = [];
CAk1 = C(x,:);
J11 = zeros(m,1);
J12 = zeros(m,1);
J22 = zeros(m,1);
for k = 1:m
	Dk = CAk1*KVL;
	dk = [Dk dk];
	DDk = [DDk zeros(k*nx,n); dk D0]; % Bk = C*A^{k-1}*K is the k-th MA coefficients matrix
	DELk = dk*dk'+D0*D0';
	Ldk = chol(DELk,'lower')\dk;
	Jk = Ldk'*Ldk;

	u = 1:k*nx;
	v = k*nx+1:k*n;
	J11(k) = maxabs(Jk(u,u));
	J12(k) = maxabs(Jk(u,v));
	J22(k) = maxabs(Jk(v,v));

	CAk1 = CAk1*A;
end
gp_qplot((1:k)',[J11 J12 J22],{'J11','J12','J22'},'set logs x');
