function [H,V,converged,relerr,niter,K,L] = cpsd_specfac(S,tol,maxi)

% Wilson's method of spectral factorization
%
% Performs a numerical inner-outer factorization of a spectral matrix, using
% Wilsons method. This implementation here is a slight modification of the
% original implemention by M. Dhamala (mdhamala@bme.ufl.edu) & G. Rangarajan
% (rangaraj@math.iisc.ernet.in), UF, Aug 3-4, 2006.
%
% modified by S K Mody (modysk@gmail.com), 22.Sepember.2016
% revised by M. Dhamala, June, 2017
% adapted for MVGC v2.0 by L. Barnett, October 2019
%
% Ref.: % The Factorization of Matricial Spectral Densities, SIAM J. Appl. Math,
% Vol. 23, No. 4, pgs 420-426 December 1972 by G T Wilson).
%
% Input:
%
% S:
%	Spectral matrix function. This should be specified as a (n x n x h)
%	array for the frequencies in the closed range [0, 0.5], divided
%	into equal intervals. The spectrum for negative frequencies is assumed
%	to be symmetric, ie:-
%		S(-f) = transpose(S(f))
%
% tol [default: 1e-7]:
%	The tolerance with which to check for convergence. Iterations stop
%	either when the number of iterations reaches a prespecified maximum
%	or when all of the following conditions are satisfied:-
%		|(L - L_prev)./L_prev| < tol
%		|(K - K_prev)./K_prev| < tol
%		|(S - K*K')./S| < tol
%	where |.| is the max norm.
%
% maxi [default: minimum of 100 and 1/sqrt(tol)]:
%
% Output:
%
% H, V:
%	H is complex array of the same size as S, and V is real symmetric
%	positive definite matrix such that for each k:
%		S(:,:,k) = H(:,:,k)*V*H(:,:,k)'
%
% K, L:
%	(n x n x h) complex array. Theoretically K is a function defined on
%	the boundary of the unit circle in the complex plane, such that:
%		S(:,:,k) = K(:,:,k)*K(:,:,k)'
%	Theoretically, K has a holomorphic extension in the complex plane to
%	all |z| < 1.). L is the upper triangular matrix that is the value of
%	K at the origin. V is related to L by:
%		V = L*L'
%
% converged:
%	Boolean value indicating whether the iteration converged to within the
%	specified tolerance.
%
% relerr:
%	The relative Cauchy error of the convergence of the spectrum or K.

	deftol  = 1e-7;
	maxmaxi = 100;

	if nargin < 2 || isempty(tol),  tol  = deftol;  end
	if nargin < 3 || isempty(maxi), maxi = min(maxmaxi,round(1/sqrt(tol))); end

	[n,~,h] = size(S);

	SS = cat(3,S,conj(S(:, :, h-1:-1:2)));
	L = L_initial__(SS);
	K = repmat(L,[1,1,h]);
	K = cat(3,K,conj(K(:,:,h-1:-1:2)));
	m = size(SS,3);

	I = eye(n);

	U = zeros(size(SS));
	for j = 1:m
		U(:,:,j) = chol(SS(:,:,j));
	end

	niter = 0;
	converged = false;
	g = zeros(n,n,m);
	while niter < maxi && ~converged
		for k = 1:m
			% Equivalent to:
			% g(:,:,k) = K(:,:,k)\SS(:,:,k)/K(:,:,k)' + I;
			W = K(:,:,k)\U(:,:,k)';
			g(:,:,k) = W*W' + I;
		end

		[gp,gp0] = PlusOperator(g);
		T = -tril(gp0,-1);
		T = T-T';

		K_prev = K;
		for k = 1:m,
			K(:,:,k) = K(:,:,k)*(gp(:,:,k) + T);
		end

		L_prev = L;
		L = L*(gp0 + T);

		% Relative cauchy error. Check on S is expensive, so check L first, then K and only then S.
		[converged relerr] = check_converged_K__(K,K_prev,L,L_prev,tol);
		if converged
			% Uncomment this next line to check for relative cauchy error in spectrum.
			[converged relerr] = check_converged_S__(SS,K,tol);
		end

		niter = niter+1;
	end

	H = zeros(n,n,h);
	for k = 1:h
		H(:,:,k) = K(:,:,k)/L;
	end

	K = K(:,:,1:h);
	V = L*L';
end

function L = L_initial__(SS)

	[n,~,m] = size(SS);

	% perform ifft to obtain gammas.
	SS = reshape(SS,[n*n,m]);
	gamma = ifft(SS.');
	gamma0 = gamma(1,:);
	gamma0 = reshape(gamma0,[n n]);

	% Remove any assymetry due to rounding error.
	% This also will zero out any imaginary values
	% on the diagonal - real diagonals are required for cholesky.
	gamma0 = real((gamma0 + gamma0')/2);
	L = chol(gamma0);

end

function [gp,gp0] = PlusOperator(g) % [ ]+operation

	[n,~,m] = size(g);
	h = ceil((m+1)/2);

	g = reshape(g, [n*n,m]);
	gammma = real(ifft(g.'));
	gammma = reshape(gammma.',[n,n,m]);

	% Take half of the zero lag
	gammma(:,:,1) = 0.5*gammma(:,:,1);
	gp0 = gammma(:,:,1);

	% Zero out negative powers.
	gammma(:,:,h+1:end) = 0;

	% Reconstitute
	gammma = reshape(gammma,[n*n,m]);
	gp = fft(gammma.');
	gp = reshape(gp.',[n,n,m]);

end

function [converged_K,relerr] = check_converged_K__(K,K_prev,L,L_prev,tol)

	[converged_K,relerr] = CheckRelErr__(L,L_prev,tol);
	if converged_K
		[converged_K RelErr2] = CheckRelErr__(K,K_prev,tol);
		relerr = max(relerr, RelErr2);
	end

end

function [converged_S,relerr] = check_converged_S__(S,K,tol)

	FX = zeros(size(K));
	for j = 1:size(K,3)
		FX(:,:,j) = K(:,:,j)*K(:,:,j)';
	end
	[converged_S,relerr] = CheckRelErr__(FX,S,tol);

end

function [ok,relerr] = CheckRelErr__(A,B,reltol)

	D = abs(B-A);
	A = abs(A);
	A(A <= 2*eps) = 1; % Minimum detectable difference between x and a value close to x is O(x)*eps.
	E = D./A;
	relerr = max(E(:));
	ok = (relerr <= reltol);

end
