function c = fracint_coeffs(d,r,arform)

% Fractional integration coefficients in AR or MA form

c = zeros(1,r+1);
c(1) = 1;
if arform % return AR form
	% Calculate the AR coefficients: a(k) = (-1)^k d(d-1)...(d-k+1)/k!, k = 0,...,r
	dd = 1+d;
else % return MA form
	% Calculate the MA coefficients: b(k) = d(d+1)...(d+k-1)/k!, k = 0,...,r
%	k = (0:r);
%	c = exp(gammaln(k+d)-gammaln(k+1)-gammaln(d));
	dd = 1-d;
end
for k = 1:r
	c(k+1) = (1-dd/k)*c(k);
end
