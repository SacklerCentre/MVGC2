function [F,pval,A,C,K,V] = tsdata_to_ss_pwcgc_permtest(X,sspf,ssmo,nperms,dclags)

% Estimate pairwise-conditional state-space GC, and calculate p-value using an
% empirical null permutation distribution. Permutations are formed by rotating
% the source channel data; rotations are offset by at least 'dclags', which
% should be set sufficiently large to ensure permuted source is decorrelated
% from target and conditioning data (see tsdata_rotate.m, decorrlags.m).
%
% Input:
%
% X        time-series data (channels x observations x trials)
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

n = size(X,1);

F = nan(n);
[A,C,K,V] = tsdata_to_ss(X,sspf,ssmo); % estimate SS model
L = chol(V,'lower'); KL  = K*L; KVK = KL*KL'; LDV = log(diag(V));
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
    [~,VR] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
    F(r,y) = log(diag(VR))-LDV(r);
end

pval = nan(n);
X0 = X;
F0 = zeros(n-1,i);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
	for i = 1:nperms
		X0(y,:,:) = tsdata_rotate(X(y,:,:),dclags);    % randomly permute (rotate) source channel time series
		[A,C,K,V]  = tsdata_to_ss(X0,sspf,ssmo);       % estimate permutation null SS model
		L = chol(V,'lower'); KL  = K*L; KVK = KL*KL';
		LDV = log(diag(V));
		[~,VR] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
		F0(:,i) = log(diag(VR))-LDV(r);
	end
	pval(r,y) = mean(F(r,y) <= F0,2);                  % p-value of F with respect to empirical permutation null distribution F0
end
