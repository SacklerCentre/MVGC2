function [FFT,FLR] = var2mvgc(A,V,x,y,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r  = [x z];              % indices of reduced variables
xr = 1:length(x);        % index of x in reduced variables

[~,VR]  = tsdata_to_var(X(r,:,:),p,regmode);  % reduced regression

FFT  = trace(VR(xr,xr))/trace(V(x,x)) - 1;    % F test statistic
FLR  = logdet(VR(xr,xr)) - logdet(V(x,x));    % likelihood-ratio test statistic
