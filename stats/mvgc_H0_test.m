function [sig,pcrit,cval] = mvgc_H0_test(pval,alpha,mhtc,tstat,nx,ny,nz,p,m,N)

% Test null hypothesis H0 of zero Granger causality

pcrit = mhtcorrect(pval,alpha,mhtc); % crtitical p-value at given significance level with multiple hypotheses correction

sig = 0+(pval < pcrit+eps); % reject H0 when sig is true (0+ converts to double)
sig(isnan(pval)) = NaN;  % ensure NaNs in the right places!

if nargout > 1 % want critical values too

	assert(nargin == 10,'For critical values, test type, problem dimensions, and number of observations and trials must be supplied');

	if strcmpi(tstat,'F')
		ftest = true;
	elseif strcmpi(tstat,'LR')
		ftest = false;
	else
		error('Unknown test statistic');
	end

	n = nx+ny+nz;
	d = p*nx*ny; % Degrees of freedom
	M = N*(m-p); % effective number of observations
	if ftest
		d2 = nx*(M-p*n)-1; % F df2
		sf = d2/d;         % F scaling factor
		cval = finv(1-pcrit,d,d2)/sf;
	else
		sf = M;            % chi^2 scaling factor = effective number of observations
		cval = chi2inv(1-pcrit,d)/sf;
	end

end
