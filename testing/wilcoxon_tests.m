%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-sided nonparametric tests - Wilcoxon (paired) and Mann-Whitney U (unpaired) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pval,esize,stat] = wilcoxon_tests(x,y,paired_samples)

% NOTE: positive effec size says x > y !!!

assert(isvector(x) && isvector(y),'Samples must be vectors');

x = x(:); % ensure column vector
y = y(:); % ensure column vector

if paired_samples % Wilcoxon signed-rank test

	n = length(x);
	assert(length(y) == n,'For paired test, samples must be the same size');

	[pval,~,stats] = signrank(x,y);

	stat = stats.signedrank; % Wilcoxon signed-rank

	maxstat = (n*(n+1))/2;

else              % Mann-Whitney U

	nx = length(x);
	ny = length(y);

	[pval,~,stats] = ranksum(x,y);

	stat = stats.ranksum - (nx*(nx+1))/2; % Mann-Whitney U

	maxstat = nx*ny;

end

esize = 2*(stat/maxstat)-1; % normalised to lie between -1 and 1
