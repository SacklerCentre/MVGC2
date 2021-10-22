function stat = var_to_gwcggc_tstat(X,V,groups,p,regmode,tstat)

% Groupwise-conditional global GC test statistics (F or likelihood ratio).
%
% Groups supplied as a cell vector of index vectors.
%
% Test statistic is calculated for each group, conditioning on all other
% variables in the system specified by the time series X.
%
% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!
%
% See stats/mvgc_pval for statistical inference. Parameters should be nx = 1,
% ny = size of group - 1, nz = number of variables - size of group.

n = size(X,1);
if isempty(V)
	[~,V]  = tsdata_to_var(X,p,regmode); % full regression
else
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match time series');
end

if strcmpi(tstat,'F')
	ftest = true;
elseif strcmpi(tstat,'LR')
	ftest = false;
else
	error('Unknown test statistic');
end

g = check_group(groups,n);

DV = diag(V);
stat = nan(g,1);
for a = 1:g
	x = groups{a};
    z = 1:n; z(x) = [];
    nx = length(x);
	VRx  = zeros(nx,1);
	for i = 1:nx
		r = [x(i) z];
		[~,VR] = tsdata_to_var(X(r,:,:),p,regmode);
		VRx(i) = VR(1,1);

	end
	if ftest
		stat(a) = sum(VRx)/sum(DV(x)) - 1;         % F-test statistic
	else
		stat(a) = sum(log(VRx)) - sum(log(DV(x))); % likelihood-ratio test statistic
	end
end
