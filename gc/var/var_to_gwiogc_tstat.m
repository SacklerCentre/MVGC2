function stat = var_to_gwiogc_tstat(X,V,groups,inout,p,regmode,tstat)

% In/out group GC test statistics (F or likelihood ratio).
%
% Groups supplied as a cell vector of index vectors.
%
% For each group, test statistic for GC is calculated either to or from the
% group and the rest of the system.
%
% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!
%
% See stats/mvgc_pval for statistical inference. Parameters should be:
% in:  nx = size of group, ny = number of variables - (size of group), nz = 0
% out: nx = number of variables - (size of group), ny = size of group, nz = 0

n = size(X,1);
if isempty(V)
	[~,V]  = tsdata_to_var(X,p,regmode); % full regression
else
	[n1,n2] = size(V);
	assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match time series');
end

if strcmpi(inout,'in')
	gcin = true;
elseif strcmpi(inout,'out')
	gcin = false;
else
	error('in/out parameter must be ''in'' or ''out''');
end

if strcmpi(tstat,'F')
	ftest = true;
elseif strcmpi(tstat,'LR')
	ftest = false;
else
	error('Unknown test statistic');
end

g = check_group(groups,n);

stat = nan(g,1);
for a = 1:g
	if gcin
		x = groups{a};
		y = 1:n; y(x) = [];
	else
		y = groups{a};
		x = 1:n; x(y) = [];
	end
	[~,VR] = tsdata_to_var(X(x,:,:),p,regmode);
	if ftest
		stat(a) = trace(VR)/trace(V(x,x)) - 1; % F-test statistic
	else
		stat(a) = logdet(VR) - logdet(V(x,x)); % likelihood-ratio test statistic
	end
end
