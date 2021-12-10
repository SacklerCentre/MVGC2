%%%%%%%%%%%%%%%%%%%%%%%%%############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Illustrates VAR GC inference between conditions for multi-trial data. The
% bootstrap samples trials with replacement.
%
% Statistical inference uses a Mann-Whitney U-test (a non-parametric unpaired
% t-test) for "statistical dominance".
%
% Note 1: this is also useful if the data is NOT multi-trial, since it may be
% artificially segmented into trials.
%
% Note 2: this demo may be adapted for, e.g., pairwise-conditional GC; in this
% case, a multiple-hypothesis correction should be applied (see 'significance.m').
%
% Note 3: statistical inference may be compromised if the VAR models have
% different model order and/or number of observations, and the total number of
% observations is small, since the bias of the GC may then be significantly
% different between conditions. A possible resolution in this case would be to
% use the likelihood-ratio (LR) test statistic rather than the single-regression
% estimator as in this routine, and to de-bias the LR estimate prior to inference
% (see 'var_to_mvgc_tstat.m' and 'mvgc_bias.m').
%
%%%%%%%%%%%%%%%%%%%%%%%%%############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters - override on command line

if ~exist('nx',     'var'), nx      = 3;           end % number of target variables
if ~exist('ny',     'var'), ny      = 5;           end % number of source variables
if ~exist('nz',     'var'), nz      = 2;           end % number of conditioning variables
if ~exist('p',      'var'), p       = [4,6];       end % model orders
if ~exist('m',      'var'), m       = [50,60];     end % numbers of observations per trial
if ~exist('N',      'var'), N       = [100,120];   end % numbers of trials
if ~exist('rho',    'var'), rho     = [0.9,0.95];  end % spectral radii
if ~exist('wvar',   'var'), wvar    = [0.9,0.7];   end % var coefficients decay weighting factors
if ~exist('rmi',    'var'), rmi     = [0.8,1.2];   end % residuals log-generalised correlations (multi-information):
if ~exist('regm',   'var'), regm    = 'OLS';       end % VAR model estimation regression mode ('OLS' or 'LWR')
if ~exist('tstat',  'var'), tstat   = 'LR';        end % GC test statistic: F or LR (likelihood ratio)
if ~exist('debias', 'var'), debias  = true;        end % Debias GC statistics? (recommended for inference)
if ~exist('alpha',  'var'), alpha   = 0.05;        end % Significance level
if ~exist('S',      'var'), S       = [1100,900];  end % bootdstrap sample sizes
if ~exist('hbins',  'var'), hbins   = 50;          end % histogram bins
if ~exist('seed',   'var'), seed    = 0;           end % random seed (0 for unseeded)
if ~exist('fignum', 'var'), fignum  = 1;           end % figure number

%%%%%%%%%%%%%%%%%%%%%%%%%############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
x = 1:nx;
y = nx+1:nx+ny;

rng_seed(seed);

Fa  = zeros(2,1);  % actual GC
Fs  = zeros(2,1);  % sample GC
Fb  = cell(2,1);   % bootstrap GCs
Fbm = zeros(2,1);  % bootstrap median
Fbd = zeros(2,1);  % bootstrap mean absolute deviation

for c = 1:2 % conditions 1 and 2

	% Random VAR model

	A = var_rand(n,p(c),rho(c),wvar(c));
	V = corr_rand(n,rmi(c));

	% Calculate actual GC y -> x

	Fa(c) = var_to_mvgc(A,V,x,y);

	% The multi-trial data

	X = var_to_tsdata(A,V,m(c),N(c));

	% Calculate sample GC y -> x

	[As,Vs] = tsdata_to_var(X,p(c),regm);
	Fs(c) = var_to_mvgc(As,Vs,x,y);

	% GC multi-trial bootstrap distribution

	fprintf('Calculating multi-trial bootstrap (condition %d) ',c);
	[Fb{c},et] = mvgc_var_multitrial_bootstrap(X,x,y,p(c),S(c),regm);
	fprintf(' %.2f seconds\n',et);

	% GC bootstrap median and median absolute deviation

	Fbm(c) = median(Fb{c});
	Fbd(c) = mad(Fb{c},1);
end

% Summary statistics (Note: you can use the median absolute deviations in Fbd
% as you would standard deviation, for error bars, etc.)

fprintf('\n------------------------------------------------\n');
fprintf('GC                   Condition 1    Condition 2\n');
fprintf('------------------------------------------------\n');
fprintf('Actual GC            : %6.4f         %6.4f\n', Fa(1),  Fa(2) );
fprintf('Sample GC            : %6.4f         %6.4f\n', Fs(1),  Fs(2) );
fprintf('Bootstrap GC median  : %6.4f         %6.4f\n', Fbm(1), Fbm(2));
fprintf('Bootstrap GC mad     : %6.4f         %6.4f\n', Fbd(1), Fbd(2));
fprintf('------------------------------------------------\n\n');

% Unpaired "t-test" (Mann-Whitney U-test) between bootstrap samples.
% Null hypothesis is that neither group stochastically dominates the
% other. A significant positive value means GC in Condition 2
% stochastically dominates (i.e. is "bigger than") GC in Condition 1.
%
% Note: if you are doing this for multiple GCs, you should use a multiple-
% hypothesis adjustment; see MVGC2 routine 'significance'.

z    = mann_whitney(Fb{1},Fb{2}); % z-score ~ N(0,1) under H0
pval = erfc(abs(z)/sqrt(2));      % p-value (2-tailed test)
sig  = pval <= alpha;             % significant (reject H0)?

if sig
	if z > 0
		sigstr = 'YES (Condition 2 > Condition 1)';
	else
		sigstr = 'YES (Condition 1 > Condition 2)';
	end
else
	sigstr = 'NO';
end

fprintf('z-score     : %6.4f\n',z);
fprintf('p-value     : %6.4f\n',pval);
fprintf('Significant : %s\n\n', sigstr);

% Plot histograms of bootstrap distributions

figure(fignum); clf;
histogram(Fb{1},hbins,'facecolor','g');
hold on
histogram(Fb{2},hbins,'facecolor','r');
hold off
title(sprintf('Bootstrap distributions\n'));
xlabel('GC (green = Condition 1, red = Condition 2)')

%%%%%%%%%%%%%%%%%%%%%%%%%############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fb,et] = mvgc_var_multitrial_bootstrap(X,x,y,p,S,regm)

	% A multi-trial bootstrap of single-regression VAR GC estimates
	%
	% X         multi-trial time-series data
	% x         target variable indices
	% y         source variable indices
	% p         VAR model order
	% S         number of samples
	% regm      regression mode ('OLS' or 'LWR')
	%
	% Fb        multi-trial bootstrap GC samples
	% et        elapsed time

	tic; % start timer
	N  = size(X,3);
	Fb = zeros(S,1); % the GC bootstrap samples
	S10 = round(S/10);
	for i = 1:S
		if rem(i,S10) == 0, fprintf('.'); end % progress indicator
		Xs      = X(:,:,randi(N,N,1));        % subsample trials with replacement
		[As,Vs] = tsdata_to_var(Xs,p,regm);   % estimated model
		Fb(i)    = var_to_mvgc(As,Vs,x,y);    % sample GC
	end
	et = toc; % elapsed time

end
