%% significance
%
% Return statistical significance and optionally critical p-values, corrected for multiple hypotheses tests
%
% <matlab:open('significance.m') code>
%
%% Syntax
%
%     [sig,pcrit] = significance(pval,alpha,correction,sym)
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
%     sym          p-values are in a symmetric matrix; only test upper triangle (default: false)
%
% _output_
%
%     sig          statistical significance (1 = reject H0, 0 = can't reject H0
%     pcrit        critical p-value
%
%% Description
%
% Returns significance and critical p-value based on p-values in |pval|,
% which may be a scalar,% vector or matrix, and significance level |alpha|.
% NaNs are ignored. The |correction| parameter specifies a multiple hypotheses
% test adjustment, and may be one of: |'None'|, |'Bonferroni'|, |'Sidak'|,
% |'FDR'| (false discovery rate, independent hypotheses or positively correlated
% hypotheses [1]) or |'FDRD'| (false discovery rate, arbitrary dependencies [2]).
%
% *_Note:_* |correction = 'None'| is not recommended for multiple
% hypotheses, so is _not_ the default! |'FDRD'| is generally a good choice.
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
% <mvgc_pval.html |mvgc_pval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [sig,pcrit] = significance(pval,alpha,correction,sym)

if nargin < 4 || isempty(sym), sym = false; end

if sym % for symmetric matrices, only use upper triangle
	[n,n1] = size(pval);
	assert(ismatrix(pval) && n1 == n,'p-values must be a square symmetric matrix');
	utidx = logical(triu(ones(n),1)); % logical indices of upper triangle
	[sig,pcrit] = significance(pval(utidx),alpha,correction,false);
	return
end

sig = NaN(size(pval)); % same shape as p-value array
fi  = isfinite(pval);  % index to finite entries (i.e., not NaN, Inf, etc.) - logical array

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

sig(fi) = 0+(pval < pcrit+eps); % reject H0 when sig is true (0+ converts to double)
                                % finites will be in same positions as they were in original p-value array

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
