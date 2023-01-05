function stat = var_to_iogc_tstat(X,V,inout,p,regmode,tstat)

% In/out GC test statistics per variable
%
% For each variable, the GC test statistic is calculated either to or
% from the variable and the rest of the system.
%
% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!
%
% See stats/mvgc_pval for statistical inference. Parameters should be:
% in:  nx = 1, ny = number of variables - 1, nz = 0
% out: nx = number of variables - 1, ny = 1, nz = 0

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

stat = nan(n,1);
if gcin
	for x = 1:n
		y = 1:n; y(x) = [];
		[~,VR] = tsdata_to_var(X(x,:,:),p,regmode);
		if ftest
			stat(i) = VR/V(x,x) - 1;               % F-test statistic
		else
			stat(i) = log(VR) - log(V(x,x));       % likelihood-ratio test statistic
		end
	end
else
	for y = 1:n
		x = 1:n; x(y) = [];
		[~,VR] = tsdata_to_var(X(x,:,:),p,regmode);
		if ftest
			stat(i) = trace(VR)/trace(V(x,x)) - 1; % F-test statistic
		else
			stat(i) = logdet(VR) - logdet(V(x,x)); % likelihood-ratio test statistic
		end
	end
end
