function F = var_to_mvgc(A,V,x,y)

% Note that at the moment the "single regression" GC estimate
% calculated here is only useful as an effect size, rather than
% hypothesis testing. In future, we intend to implement a
% hypothesis test for this statistic; see
%
%     A. J. Gutknecht and L. Barnett, Sampling distribution for
%     single-regression Granger causality estimators,
%     arXiv 1911.09625 [math.ST], 2019.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match VAR coefficients matrix');

x = x(:)'; % vectorise target variable indices
y = y(:)'; % vectorise source variable indices

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'Some x indices out of range');
assert(all(y >=1 & y <= n),'Some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of conditioning variables (i.e., all other variablesin model)
r = [x z];               % indices of variables in reduced model (omit source variables)

xr = 1:length(x); % indices of target in reduced model

F = NaN;
[~,VR,rep] = var2riss(A,V,y,r);
if sserror(rep), return; end % check DARE report, bail out on error
F = logdet(VR(xr,xr)) - logdet(V(x,x));
