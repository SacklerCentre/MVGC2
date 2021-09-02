%% significance
%
% Statistical significance adjusted for multiple hypotheses
%
% <matlab:open('significance.m') code>
%
%% Syntax
%
%     sig = significance(pval,alpha,correction)
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
%     sig          significances (0 or 1)
%     pcrit        critical p-value
%
%% Description
%
% Returns significance (0 or 1) of statistics based on p-values in |pval|,
% which may be a scalar, vector or matrix. NaNs are ignored. The
% |correction| parameter specifies a multiple hypotheses test adjustment,
% and may be one of: |'None'|, |'Bonferroni'|, |'Sidak'|,
% |'Holm'|, |'FDR'| (false discovery rate, independent or
% positively correlated hypotheses [1]) or |'FDRD'| (false discovery rate,
% arbitrary dependencies [2]).
%
% *_Note:_* |correction = 'None'| is not recommended for multiple
% hypotheses, so is _not_ the default! |'FDR'| is generally a good choice.
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
% <mvgc_pval.html |mvgc_pval|> |
% <empirical_pval.html |empirical_pval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [sig,pcrit] = significance(pval,alpha,correction)

sig = NaN(size(pval)); % same shape as p-value array
fi  = isfinite(pval);  % index to finite entries (i.e., not NaN, Inf, etc.) - logical array
p   = pval(fi);        % vectorise the finite p-values
n   = numel(p);        % number of p-values being tested

switch upper(correction)

    case 'NONE';

		pcrit = alpha;
		signn = (p < pcrit);

    case 'BONFERRONI' % assumes independence of test statistics

        pcrit = alpha/n;
		signn = (p < pcrit);

    case 'SIDAK' % assumes independence of test statistics

        pcrit = 1-realpow(1-alpha,1/n);
		signn = (p < pcrit);

    case 'HOLM' % v2.0 DOESN'T assume independence of test statistics!

		pcrit = NaN; % temporary (not useful)
        signn = false(1,n);
        [psorted,sortidx] = sort(p);
        for v=1:n
            if psorted(v) < alpha/(n-v+1)
                signn(sortidx(v)) = true;
            else
                break; % remaining null hypotheses accepted
            end
        end

    case 'FDR'   % assumes independence (or positive correlation) of test statistics (more powerful)

        pcrit = fdr_bh(p,alpha,true);
		signn = (p < pcrit);

    case 'FDRD' %  possible dependencies - no correlation assumptions

        pcrit = fdr_bh(p,alpha,false);
		signn = (p < pcrit);

    otherwise; error('unknown correction method');
end

sig(fi) = signn; % finites will be in same positions as they were in original p-value array

function pcrit = fdr_bh(p,q,pdep)

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

s = size(p);
if length(s) > 2 || s(1) > 1
    psorted = sort(reshape(p,1,prod(s)));
else % p-values are already a row vector
    psorted = sort(p);
end
m = length(psorted); % number of tests

if pdep % BH procedure for independence or positive dependence
    thresh = (1:m)*q/m;
else    % BH procedure for any dependency structure
    thresh = (1:m)*q/(m*sum(1./(1:m)));
end

rej = psorted <= thresh;
max_id = find(rej,1,'last'); % find greatest significant pvalue
if isempty(max_id),
    pcrit = 0;
else
    pcrit = psorted(max_id);
end
