function F = var_to_pwcgc(A,V)

% Note that at the moment the "single regression" GC estimates
% calculated here are only useful as an effect size, rather than
% hypothesis testing. In future, we intend to implement a
% hypothesis test for this statistic; see
%
%     A. J. Gutknecht and L. Barnett, Sampling distribution for
%     single-regression Granger causality estimators,
%     arXiv 1911.09625 [math.ST], 2019.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

DV = diag(V);
LDV = log(DV);

F = nan(n);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
	[~,VR,rep] = var2riss(A,V,y,r);
    if sserror(rep,y), continue; end % check DARE report, bail out on error
    F(r,y) = log(diag(VR))-LDV(r);
end
