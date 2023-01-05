function [F,et] = mvgc_var_bootstrap(X,x,y,p,S,regm)

% A multi-trial bootstrap of single-regression VAR GC estimates
%
% X         multi-trial time-series data
% x         target variable indices
% y         source variable indices
% p         VAR model order
% S         number of samples
% regm      regression mode ('OLS' or 'LWR')
%
% F         multi-trial bootstrap GC test statistic sample
% et        elapsed time

tic;
N = size(X,3);
s = randi(N,N,S); % subsample trials with replacement
F = zeros(S,1);
S10 = round(S/10);
for i = 1:S
	if rem(i,S10) == 0, fprintf('.'); end % progress indicator
	Xs      = X(:,:,s(:,i));              % select i-th bootstrap sample
	[As,Vs] = tsdata_to_var(Xs,p,regm);   % estimated model
	F(i)    = var_to_mvgc(As,Vs,x,y);     % sample GC test statistic
end
et = toc; % elapsed time
