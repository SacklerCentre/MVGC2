function [Y,X,mtrunc] = varfima_to_tsdata(A,B,fiparms,X,m,N,mtrunc,decfac)

% [Y,X,mtrunc] = varfima_to_tsdata(A,B,fiparms,V,m,N,mtrunc,decfac)
%
% or
%
% [Y,X,mtrunc] = varfima_to_tsdata(A,B,fiparms,X,[],[],mtrunc,decfac)

have_var  = ~isempty(A);
have_vma  = ~isempty(B);
fint      = ~isempty(fiparms);

if fint
	assert(isvector(fiparms) && length(fiparms) == 3,'fractional integration parameters must be a 2-vector [d,r]');
	d = fiparms(1);
	r = fiparms(2);
	arform = fiparms(3);
end

if nargin < 7 || isempty(mtrunc) % automatic calculation - transients decay with rate given by VAR spectral radius
    if nargin < 8 || isempty(decfac), decfac = 1; end
	if have_var
		if have_vma
			rho = max(specnorm(A),specnorm(-B));
		else
			rho = specnorm(A);
		end
	else
		if have_vma
			rho = specnorm(-B);
		else
			rho = NaN;
		end
	end
	if isnan(rho)
		mtrunc = 0;
	else
		mtrunc = ceil(decfac*(-log(eps))/(-log(rho)));
		if fint
			fprintf(2,'WARNING: automatic truncation doesn''t really work for fractional integration!\n');
		end
	end
else
    assert(isscalar(mtrunc) && isint(mtrunc) && mtrunc >= 0,'truncation parameter must be a non-negative integer');
end

[n,m1,N1] = size(X);
if nargin < 5 || isempty(m) % X is input
	assert(nargin < 6 || isempty(N),'number of trials specified by input!');
	m = m1;
	N = N1;
else          % X is a covariance matrix - generate Gaussian residuals as input
	if nargin < 6 || isempty(N), N = 1; end % single trial
	assert(ismatrix(X) && m1 == n,'covariance matrix not square');
	[L,cholp] = chol(X,'lower');
	assert(cholp == 0,'covariance matrix not positive-definite');
	m = m+mtrunc;
	X = zeros(n,m,N);
	for k = 1:N
		X(:,:,k) = L*randn(n,m); % Gaussian residuals
	end
end

assert(mtrunc < m,'too much truncation');

if N > 1 % multi-trial
    Y = zeros(n,m,N);
  	if have_var
		if have_vma
			if fint
				for k = 1:N
					Y(:,:,k) = gendiff(genvar(A,genvma(B,X(:,:,k))),fiparms);
				end
			else
				for k = 1:N
					Y(:,:,k) = genvar(A,genvma(B,X(:,:,k)));
				end
			end
		else
			if fint
				for k = 1:N
					Y(:,:,k) = gendiff(genvar(A,X(:,:,k)),fiparms);
				end
			else
				for k = 1:N
					Y(:,:,k) = genvar(A,X(:,:,k));
				end
			end
		end
	else
		if have_vma
			if fint
				for k = 1:N
					Y(:,:,k) = gendiff(genvma(B,X(:,:,k)),fiparms);
				end
			else
				for k = 1:N
					Y(:,:,k) = genvma(B,X(:,:,k));
				end
			end
		else
			if fint
				for k = 1:N
					Y(:,:,k) = gendiff(X(:,:,k),fiparms);
				end
			else
				Y = X;
			end
		end
	end
else
	if have_var
		if have_vma
			if fint
				Y = gendiff(genvar(A,genvma(B,X)),fiparms);
			else
				Y = genvar(A,genvma(B,X));
			end
		else
			if fint
				Y = gendiff(genvar(A,X),fiparms);
			else
				Y = genvar(A,X);
			end
		end
	else
		if have_vma
			if fint
				Y = gendiff(genvma(B,X),fiparms);
			else
				Y = genvma(B,X);
			end
		else
			if fint
				Y = gendiff(X,fiparms);
			else
				Y = X;
			end
		end
	end
end

if mtrunc > 0
	Y = Y(:,mtrunc+1:m,:);
	if nargout > 1
		X = X(:,mtrunc+1:m,:);
	end
end
