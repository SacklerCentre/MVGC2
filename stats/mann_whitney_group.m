function [z,U] = mann_whitney_group(x1,x2)

% "Groupwise" Mann-Whitney U test - a nonparametric unpaired "t-test" for
% stochastic dominance between independent paired groups of samples; that
% is, comparisons are only between samples belonging to the same group-pair.
% Positive effect means 2nd group stochastically dominates 1st group.
%
% Typical usage is for a cross-subject analysis, where each "group" corresponds
% to a subject.
%
% x1    first group of sets of values (cell-array of vectors)
% x2    second group of sets of values (cell-array of vectors, same number of cells as x1)
%
% z     z-score (asymptotically standard normal)
% U     the Mann-Whitney U-statistic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To calculate p-value(s)
%
% p = erfc(abs(z)/sqrt(2)); % two-tailed
% p = erfc(+z/sqrt(2))/2;   % right-tail
% p = erfc(-z/sqrt(2))/2;   % left-tail
%
% This works also for multiple hypothesis testing, in which case z may be a vector).
% To calculate significances and critical z-score at significance level alpha with
% multiple-hypothesis adjustment mhtc:
%
% [sigs,pcrit] = significance(pvals,alpha,mhtc);
%
% then
%
% zcrit = (+/-)sqrt(2)*erfcinv(pcrit); % two-tailed (+/-)
% zcrit = +sqrt(2)*erfcinv(2*pcrit);   % right-tail
% zcrit = -sqrt(2)*erfcinv(2*pcrit);   % left-tail
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(iscell(x1) && iscell(x2) && isvector(x1) && isvector(x2) && length(x1) == length(x2),'Paired-group data must be matching cell vectors of samples');

% Calculate Mann-Whitney U for each group

N  = length(x1); % number of groups
u  = zeros(N,1);
n1 = zeros(N,1);
n2 = zeros(N,1);
for g = 1:N % for each group
	n1(g) = numel(x1{g});              % size of 1st group
	n2(g) = numel(x2{g});              % size of 2nd group
	u(g)  = nnz(x2{g}(:) > x1{g}(:)'); % n1{g}*n2{g} comparisons: count how many times x2 > x1 (ignore ties; this is floating-point!)
end

% Aggregate Mann-Whitney U across groups (groups assumed independent!)

U = sum(u);

% Calculate Mann-Whitney aggregated U asymptotic null mean, variance and z-score

m = sum((n1.*n2)/2);       % aggregated U asymptotic mean under null
v = sum((m.*(n1+n2+1))/6); % aggregated U asymptotic variance under null
z = (U-m)/sqrt(v);         % z-score ~ N(0,1) under H0
