function [F,pval,A,V] = tsdata_to_var_mvgc_permtest(X,x,y,varmo,regmode,nperms,dclags)

% Estimate conditional VAR/SS GC, and calculate p-value using an empirical
% null permutation distribution. Permutations are formed by rotating the source
% channel data; rotations are offset by at least 'dclags', which should be set
% sufficiently large to ensure permuted source is decorrelated from target and
% conditioning data (see tsdata_rotate.m, decorrlags.m).
%
% Input:
%
% X        time-series data (channels x observations x trials)
% x        target channel indices
% y        source channel indices (remaining channels are conditioned out)
% varmo    VAR model order
% regmode  VAR model regression mode: 'OLS' or 'LWR'
% nperms   number of permutations
% dclags   decorrelation lags for permutations (see decorrlags.m for estimation)
%
% Output:
%
% F        conditional GC estimate
% pval     p-value for estimate using permutation empirical null distribution
% A,V      estimated VAR parameters

[A,V] = tsdata_to_var(X,varmo,regmode); % estimate VAR model
F = var_to_mvgc(A,V,x,y);               % calculate GC

X0 = X;
F0 = zeros(nperms,1);
for i = 1:nperms
	X0(y,:,:) = tsdata_rotate(X(y,:,:),dclags); % randomly permute (rotate) source channel time series
	[A0,V0] = tsdata_to_var(X0,varmo,regmode);  % estimate null VAR model
	F0(i)   = var_to_mvgc(A0,V0,x,y);           % calculate GC
end
pval = mean(F <= F0); % p-value of F with respect to empirical permutation null distribution F0
