function stat = var_to_gwcgc_tstat(X,V,groups,p,regmode,tstat)

% Groupwise-conditional GC test statistics (F or likelihood ratio).
%
% Groups supplied as a cell vector of index vectors.
%
% Test statistic is calculated for each pair of groups (first index is target
% group, second is source group), conditioning on all other variables in the
% system specified by the time series X.
%
% NOTE: If full-regression residuals covariance matrix V is supplied, it must
% have been obtained using 'tsdata_to_var' with the SAME 'p' and  'regmode' !!!
%
% See stats/mvgc_pval for statistical inference. Parameters should be:
% nx = size of target group, ny = size of source group,
% nz = number of variables - (size of target group + size of source group).

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

stat = nan(g);
for b = 1:g
	y = groups{b};
    r = 1:n; r(y) = []; % omit group b
	[~,VR] = tsdata_to_var(X(r,:,:),p,regmode);
    for a = 1:g
        if a == b, continue; end
        x = groups{a};
        xr = findin(x,r); % indices of group{a} in r
		if ftest
			stat(a,b) = trace(VR(xr,xr))/trace(V(x,x)) - 1; % F-test statistic
		else
			stat(a,b) = logdet(VR(xr,xr))-logdet(V(x,x));   % likelihood-ratio test statistic
		end
    end
end
