function [z,U] = mann_whitney(x1,x2)

% Mann-Whitney U test - a nonparametric unpaired "t-test" for
% stochastic dominance. Positive effect means 2nd sample
% stochastically dominates 1st sample.
%
% x1    first sample  (vector)
% x2    second sample (vector)
%
% z     z-score (asymptotically standard normal)
% U     the Mann-Whitney U-statistic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To calculate p-value(s) for a two-tailed significance test (this works also for
% multiple hypothesis testing, in which case z may be a vector):
%
% pvals = erfc(abs(z)/sqrt(2));
%
% To calculate significances and critical z-score (two-tailed test) at significance
% level alpha with multiple-hypothesis adjustment mhtc:
%
% [sigs,pcrit] = significance(pvals,alpha,mhtc);
% zcrit = sqrt(2)*erfcinv(pcrit);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Mann-Whitney U

n1 = length(x1);          % size of 1st sample
n2 = length(x2);          % size of 2nd sample
U = nnz(x2(:) > x1(:)');  % n1*n2 comparisons: count how many times x2 > x1 (ignore ties; this is floating-point!)

% Calculate Mann-Whitney U asymptotic null mean, variance and z-score

m = (n1*n2)/2;            % u asymptotic mean under null
v = (m*(n1+n2+1))/6;      % u asymptotic variance under null
z = (U-m)/sqrt(v);        % z-score ~ N(0,1) under H0
