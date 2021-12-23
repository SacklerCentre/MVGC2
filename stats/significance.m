%% significance
%
% Return statistical significance and optionally critical p-values, corrected for multiple hypotheses tests
%
% <matlab:open('significance.m') code>
%
%% Syntax
%
%     [sig,pcrit] = significance(pval,alpha,correction,flag)
%     sig         = significance(nhyp,alpha,correction,'nopv')
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     pval             p-values
%     nhyp             number of hypotheses
%     alpha            significance level
%     correction       multiple hypotheses correction (see Description)
%     flag ==  empty   p-values for all finite entries taken into account
%     flag == 'symm'   only finite p-values on and above the diagonal taken into account
%     flag == 'nopv'   no p-values: first argument is number of hypotheses, only return critical p-value
%
% _output_
%
%     sig              statistical significance (1 = reject H0, 0 = can't reject H0
%     pcrit            critical p-value
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
% If the 'nopv' flag is supplied, the first argument is the number of hypothesis,
% and only critical values are returned; this may not work for some multiple
% hypothesis correction methods.
%
% If the 'symm' flag is supplied, a symmetric matrix of p-values is assumed, and
% only finite p-values on and above the diagonal are taken into account. This is
% useful, e.g., for naturally symmetric statistics such as correlation matrices.
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

function [oarg1,oarg2] = significance(iarg1,alpha,correction,flag)

nopv = false;
symm = false;
if nargin > 3 && ~isempty(flag)
	switch lower(flag)
		case 'nopv', nopv = true;
		case 'symm', symm = true;
		otherwise, error('Unknown flag');
	end
end

if nopv % some methods don't require actual p-values, just the number of hypotheses
	assert(nargout < 2,'"No p-values" option only returns one argument (critical p-value)!');
	nhyp = iarg1;  % first input argument is the number of hypotheses
else
	pval = iarg1; % first input argument is array of p-values
	if symm
		assert(ismatrix(pval),'p-values must be a matrix (assumed symmetric) if called with ''symm'' flag');
		pval(logical(tril(ones(size(pval,1)),-1))) = NaN; % NaNs below the diagonal
	end
	oarg1 = NaN(size(pval));      % first return argument is significances - same shape as p-value array
	fi    = isfinite(pval);       % index to finite entries (i.e., not NaN, Inf, etc.) - logical array
	pval  = pval(isfinite(pval)); % vectorise the finite p-values (i.e., not NaN, Inf, etc.)
	nhyp  = numel(pval);          % number of hypotheses
end

switch upper(correction)
    case 'NONE'
		pcrit = alpha;
    case 'BONFERRONI' % assumes independence of test statistics
		pcrit = alpha/nhyp;
    case 'SIDAK'      % assumes independence of test statistics
		pcrit = 1-realpow(1-alpha,1/nhyp);
    case 'FDR'        % assumes independence (or positive correlation) of test statistics (more powerful)
		assert(~nopv,'Need actual p-values for this method');
		pcrit = fdr_bh(pval,alpha,true);
    case 'FDRD'       %  possible dependencies - no correlation assumptions
		assert(~nopv,'Need actual p-values for this method');
		pcrit = fdr_bh(pval,alpha,false);
    otherwise, error('Unknown correction method');
end

if nopv
	oarg1 = pcrit;                    % first (and only) return argument is critical value
else
	oarg1(fi) = 0+(pval < pcrit+eps); % first return argument is significance; reject H0 when true (0+ converts to double), finites will be in same positions as they were in original p-value array
	oarg2 = pcrit;                    % second return argument is critical value
	if symm
		oarg1 = symmetrise(oarg1);    % copy significances from upper to lower triangle of matrix
	end
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

maxidx = find(psorted <= thresh,1,'last'); % find greatest significant pvalue
if isempty(maxidx),
    pcrit = 0;
else
    pcrit = psorted(maxidx);
end
