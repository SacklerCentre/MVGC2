function [F,pval,A,C,K,V] = tsdata_to_ss_mvgc_permtest(X,x,y,sspf,ssmo,nperms,dclags)

% Estimate conditional state-space GC, and calculate p-value using an empirical
% null permutation distribution. Permutations are formed by rotating the source
% channel data; rotations are offset by at least 'dclags', which should be set
% sufficiently large to ensure permuted source is decorrelated from target and
% conditioning data (see tsdata_permute.m, decorrlags.m).
%
% Input:
%
% X        time-series data (channels x observations x trials)
% x        target channel indices
% y        source channel indices (remaining channels are conditioned out)
% sspf     past/future horizons for CCA state-space-subspace algorithm
% ssmo     state-space model order
% nperms   number of permutations
% dclags   decorrelation lags for permutations (see decorrlags.m for estimation)
%
% Output:
%
% F        conditional GC estimate
% pval     p-value for estimate using permutation empirical null distribution
% A,C,K,V  estimated ISS parameters

[A,C,K,V] = tsdata_to_ss(X,sspf,ssmo); % estimate SS model
F = ss_to_mvgc(A,C,K,V,x,y);           % calculate GC

X0 = X;
F0 = zeros(nperms,1);
for i = 1:nperms
	X0(y,:,:) = tsdata_permute(X(y,:,:),dclags); % randomly permute (rotate) source channel time series
	[A,C,K,V] = tsdata_to_ss(X0,sspf,ssmo);      % estimate permutation null SS model
	F0(i)     = ss_to_mvgc(A,C,K,V,x,y);         % sample GC estimate from permutation null distribution
end
pval = mean(F <= F0); % p-value of F with respect to empirical permutation null distribution F0
