function [F,pval,A,V] = tsdata_to_var_pwcgc_permtest(X,varmo,regmode,nperms,dclags)

% Estimate pairwise-conditional VAR/SS GC, and calculate p-value using an
% empirical null permutation distribution. Permutations are formed by rotating
% the source channel data; rotations are offset by at least 'dclags', which
% should be set sufficiently large to ensure permuted source is decorrelated
% from target and conditioning data (see tsdata_permute.m, decorrlags.m).
%
% Input:
%
% X        time-series data (channels x observations x trials)
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

n = size(X,1);

F = nan(n);
[A,V] = tsdata_to_var(X,varmo,regmode); % estimate VAR model
LDV = log(diag(V));
for y = 1:n
	r = [1:y-1 y+1:n]; % omit y
	[~,VR] = vardare(A,V,y,r); % "reduced" innovations covariance
	F(r,y) = log(diag(VR))-LDV(r);
end

pval = nan(n);
X0 = X;
F0 = zeros(n-1,i);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
	for i = 1:nperms
		X0(y,:,:) = tsdata_permute(X(y,:,:),dclags); % randomly permute (rotate) source channel time series
		[A,V] = tsdata_to_var(X,varmo,regmode);      % estimate VAR model
		LDV = log(diag(V));
		[~,VR] = vardare(A,V,y,r);                   % "reduced" innovations covariance
		F0(:,i) = log(diag(VR))-LDV(r);
	end
	pval(r,y) = mean(F(r,y) <= F0,2);                % p-value of F with respect to empirical permutation null distribution F0
end
