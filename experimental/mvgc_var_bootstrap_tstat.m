function [F,et] = mvgc_var_bootstrap_tstat(X,x,y,p,S,regm,tstat)

% A multi-trial bootstrap of F or likelihood ratio VAR GC test statistics
%
% NOTE: the advantage of this over using the single-regression estimator
% as in mvgc_var_bootstrap_tstat.m is that the test statistics may be
% approximately de-biased; see mvgc_bias.m.
%
% X         multi-trial time-series data
% x         target variable indices
% y         source variable indices
% p         VAR model order
% S         number of samples
% regm      regression mode ('OLS' or 'LWR')
% tstat     test statistic: F or likelihood ratio (chi^2)
%
% F         multi-trial bootstrap GC sample
% et        elapsed time

tic;
N = size(X,3);
s = randi(N,N,S); % subsample trials with replacement
F = zeros(S,1);
S10 = round(S/10);
for i = 1:S
	if rem(i,S10) == 0, fprintf('.'); end             % progress indicator
	Xs   = X(:,:,s(:,i));                             % select i-th bootstrap sample
	F(i) = var_to_mvgc_tstat(Xs,[],x,y,p,regm,tstat); % GC test statistic (F or likelihood-ratio)
end
et = toc; % elapsed time
