function [z,U] = mann_whitney_group(x1,x2)

% "Groupwise" Mann-Whitney U test - a nonparametric unpaired "t-test"
% for stochastic dominance between paired groups of samples; that is,
% comparison is only between samples belonging to the same group.
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
% To calculate p-value(s) for a two-tailed significance test (this works also for
% multiple hypothesis testing, in which case z may be a vector):
%
% pvals = erfc(abs(z)/sqrt(2));
%
% To calculate significances and critical z-score at significance level alpha
% with multiple-hypothesis adjustment mhtc:
%
% [sigs,pcrit] = significance(pvals,alpha,mhtc);
% zcrit = sqrt(2)*erfcinv(pcrit);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Mann-Whitney U-statistics for each group

assert(iscell(x1) && iscell(x2) && isvector(x1) && isvector(x2) && length(x1) == length(x2),'Data must be matching cell vectors');

% Mann-Whitney U test statistics

N  = length(x1); % number of groups
u  = zeros(N,1);
n1 = zeros(N,1);
n2 = zeros(N,1);
for g = 1:N % for each group
	n1(g) = numel(x1{g});              % size of 1st group
	n2(g) = numel(x2{g});              % size of 2nd group
	u(g)  = nnz(x2{g}(:) > x1{g}(:)'); % n1{g}*n2{g} comparisons: count how many times x2 > x1 (ignore ties; this is floating-point!)
end

% Calculate Mann-Whitney U theoretical null means, variances and z-score under normal approximation

m = (n1.*n2)/2;               % u theoretical means under null
v = (m.*(n1+n2+1))/6;         % u theoretical variances under null
M = sum(m);                   % U mean under null
V = sum(v);                   % U variance under null (groups assumed independent!)
S = sqrt(V);                  % U standard deviation under null
U = sum(u);                   % U aggregated across groups (groups assumed independent!)
z = (U-M)/S;                  % z-score ~ N(0,1) under H0
