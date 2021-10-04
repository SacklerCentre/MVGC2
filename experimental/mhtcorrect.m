%% mhtcorrect
%
% Return critical p-value corrected for multiple hypotheses tests
%
% <matlab:open('significance.m') code>
%
%% Syntax
%
%     pcrit = mhtcorrect(pval,alpha,correction)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     pval         p-values
%     alpha        significance level
%     correction   multiple hypotheses correction (see Description)
%
% _output_
%
%     pcrit        critical p-value
%
%% Description
%
% Returns critical p-value based on p-values in |pval|, which may be a scalar,
% vector or matrix, and significance level |alpha|. NaNs are ignored. The
% |correction| parameter specifies a multiple hypotheses test adjustment,
% and may be one of: |'None'|, |'Bonferroni'|, |'Sidak'|, |'FDR'| (false
% discovery rate, independent hypotheses or positively correlated hypotheses [1])
% or |'FDRD'| (false discovery rate, arbitrary dependencies [2]).
%
% *_Note:_* |correction = 'None'| is not recommended for multiple
% hypotheses, so is _not_ the default! |'FDRD'| is generally a good choice.
%
% Given multiple p-values, multiple hypothesis tests may be performed as
%
%     psig = mhtcorrect(pval,alpha,correction);
%     h = pval <= pcrit;
%
%% References
%
% [1] Y. Benjamini and Y. Hochberg, "Controlling the
% false discovery rate: a practical and powerful approach to multiple
% testing", _J. Royal Stat. Soc. B_, 57(1), 1995.
%
% [2] Y. Benjamini and D. Yekutieli, "The control of the false discovery
% rate in multiple testing under dependency", _Ann. Stat_, 29(4), 2001.
%
%% See also
%
% <var_to_mvgc_stats.html |var_to_mvgc_stats|> |
% <empirical_pval.html |empirical_pval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function pcrit = mhtcorrect(pval,alpha,correction)

% some methods don't require actual p-values

no_pvals = (isscalar(pval) && pval < 0); % -pval is now just the number of tests

if no_pvals
	n = -pvals;
else
	pval = pval(isfinite(pval)); % vectorise the finite p-values (i.e., not NaN, Inf, etc.)
	n = numel(pval);             % number of p-values being tested
end

switch upper(correction)

    case 'NONE';

		pcrit = alpha;

    case 'BONFERRONI' % assumes independence of test statistics

		pcrit = alpha/n;

    case 'SIDAK' % assumes independence of test statistics

		pcrit = 1-realpow(1-alpha,1/n);

    case 'FDR'   % assumes independence (or positive correlation) of test statistics (more powerful)

		assert(~no_pvals,'Need actual p-values for this method');
		pcrit = fdr_bh(pval,alpha,true);

    case 'FDRD' %  possible dependencies - no correlation assumptions

		assert(~no_pvals,'Need actual p-values for this method');
		pcrit = fdr_bh(pval,alpha,false);

    otherwise
		error('Unknown correction method');
end

function pcrit = fdr_bh(pvals,q,pdep)

% Executes the Benjamini & Hochberg (1995) procedure for
% controlling the false discovery rate (FDR) of a family of
% hypothesis tests. FDR is the expected proportion of rejected
% hypotheses that are mistakenly rejected (i.e., the null
% hypothesis is actually true for those tests). FDR is a
% somewhat less conservative/more powerful method for correcting
% for multiple comparisons than methods like Bonferroni
% correction that provide strong control of the family-wise
% error rate (i.e., the probability that one or more null
% hypotheses are mistakenly rejected).
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010

psorted = sort(pvals(:)');

m = length(psorted); %number of tests
if pdep % BH procedure for independence or positive dependence
	thresh = (1:m)*q/m;
else    % BH procedure for any dependency structure
	thresh = (1:m)*q/(m*sum(1./(1:m)));
end

maxidx = find(psorted <= thresh,1,'last'); %find greatest significant pvalue
if isempty(maxidx),
    pcrit = 0;
else
    pcrit = psorted(maxidx);
end
