r   = 14;
nx  = 5;
ny  = 11;
rho = 0.96;
g   = 1;

kmax   = 1000;

tol = eps;
%tol = sqrt(eps);

%rng_seed(867817212);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny;
x = 1:nx;
y = nx+1:n;

[A,C,K] = iss_rand(n,r,rho);

[V,C0] = corr_rand(n,g);

ssinfo = ss_info(A,C,K,V);

[F1,k1,dF1] = ifgc(A,C,K,V,x,y,kmax,tol);

Mbavail = memavail;
Mbmax   = (((kmax+1)^2)*nx*n)/1e6;
Mbused  = (((k1+1)^2)*nx*n)/1e6;


fprintf('total performance\n');
fprintf('-----------------\n');
fprintf('iterations = %d',k1);
if dF1 > tol
	fprintf(2,' *** failed to converge ***\n');
else
	fprintf('\n');
end
fprintf('delta(F)   = %.4g (tol = %.4g)\n',dF1,tol);
fprintf('memory (Mb): available ~ %d, max ~ %d, used ~ %d\n\n',floor(Mbavail),round(Mbmax),round(Mbused));

[F2,k2,dF2] = fhgc(A,C,K,V,x,y,kmax,tol);
F2 = cumsum(F2);

fprintf('cumulative performance\n');
fprintf('----------------------\n');
fprintf('iterations = %d',k2);
if dF2 > tol
	fprintf(2,' *** failed to converge ***\n');
else
	fprintf('\n');
end
fprintf('delta(F)   = %.4g (tol = %.4g)\n\n',dF2,tol);

k = max(k1,k2);
F1 = [F1;nan(k-k1,1)];
F2 = [F2;nan(k-k2,1)];
gp_qplot((1:k)',[F1 F2],{'total','cumulative'},'set key bottom right');

gp_qplot((1:k)',[[NaN;diff(F1)] [NaN;diff(F2)]],{'total','cumulative'},'set key bottom right');
