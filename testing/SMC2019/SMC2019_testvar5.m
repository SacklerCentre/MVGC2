%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rho',  'var'), rho   = 0.9; end
if ~exist('w',    'var'), w     = 0.2; end
if ~exist('g',    'var'), g     = 0.3; end
if ~exist('c',    'var'), c     = 0.5; end
if ~exist('vseed','var'), vseed = 0;   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

con = [...
	1 2 11  5; ...
	2 1  5  3; ...
	3 1  8 -6; ...
	3 5  4  3; ...
	4 3 20 -17 ...
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncons = size(con,1);
n = max(max(con(:,1:2)));
p = max(con(:,3));

nmh = p*n*(n-1);  % number of hypotheses to test

A = zeros(n,n,p);
if c == 0
	A(:,:,1) = diag([0.3558 0.3558 0.3558 0.3558 0.3558]);
	A(1,2,11) =  0.2214;
	A(2,1, 5) =  0.3060;
	A(3,1, 8) = -0.4033;
	A(3,5, 4) =  0.3516;
	A(4,3,20) = -0.2154;
	specnorm(A)
else
	A(:,:,1) = eye(n);
	for i = 1:ncons
		A(con(i,1),con(i,2),con(i,3)) = c*con(i,4);
	end
	A = specnorm(exp(-w*sqrt(p))*A,rho);
end

for i = 1:ncons
	con(i,5) = A(con(i,1),con(i,2),con(i,3));
end

rstate = rng_seed(vseed);
V = corr_rand(n,g);
rng_restore(rstate);
