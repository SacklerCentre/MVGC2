function [Y,o] = tsdata_rotate(X,mino)

% Left-rotate time series X randomly by at least mino observations.
% All variables in a given trial are rotated by the same offset, but
% different trials use different offsets.
%
% For GC permutation testing, ensure offset large enough that
% autocorrelation has decayed below statistical significance.
% Autocorrelation decays roughly exponentially, with decay rate
% given by spectral radius; see function 'decorrlags'.

[nvars,nobs,ntrials] = size(X);

assert(2*mino <= nobs && mino >= 1,'Minimum offset out of range');

a = nobs-2*mino+1;
b = mino-1;
Y = zeros(nvars,nobs,ntrials);
for k = 1:ntrials
	o = randi(a)+b;  % uniform in range mino:(nobs-mino)
	Y(:,:,k) = [X(:,o+1:nobs,k) X(:,1:o,k)]; % left-rotate trial
end
