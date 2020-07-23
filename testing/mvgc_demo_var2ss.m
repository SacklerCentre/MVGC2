%% MVGC demo
%
% Demonstrates typical usage of the MVGC toolbox on generated VAR data for a
% 5-node network with known causal structure (see <var5_test.html |var5_test|>).
% Estimates a VAR model and calculates time- and frequency-domain
% pairwise-conditional Granger causalities (also known as the "causal graph").
% Also calculates Seth's causal density measure [2].
%
% This script is a good starting point for learning the MVGC approach to
% Granger-causal estimation and statistical inference. It may serve as a useful
% template for your own code. The computational approach demonstrated here will
% make a lot more sense alongside the reference document> [1], which we
% _strongly recommend_ you consult, particularly Section 3 on design principles
% of the toolbox. You might also like to refer to the <mvgc_schema.html schema>
% of MVGC computational pathways - <mvgc_schema.html#3 algorithms> |A<n>| in
% this demo refer to the algorithm labels listed there - and the
% <mvgchelp.html#4 Common variable names and data structures> section of the
% Help documentation.
%
% *_FAQ:_* _Why do the spectral causalities look so smooth?_
%
% This is because spectral quantities are calculated from the estimated VAR,
% rather than sampled directly. This is in accordance with the MVGC design
% principle that all causal estimates be based on the <mvgc_demo.html#6
% estimated VAR model> for your data, and guarantees that spectral causalities
% <mvgc_demo.html#10 integrate correctly> to time-domain causality as theory
% requires. See [1] for details.
%
% *_Note_*: Do _not_ pre-filter your data prior to GC estimation, _except_
% possibly to improve stationarity (e.g notch-filtering to eliminate line noise
% or high-pass filtering to suppress low-frequency transients). Pre-filtering
% (of stationary data) may seriously degrade Granger-causal inference! If you
% want (time-domain) GC over a limited frequency range, rather calculate
% "band-limited" GC; to do this, calculate frequency-domain GCs over the full
% frequency range, then integrate over the desired frequency band [3]; see
% <smvgc_to_mvgc.html |smvgc_to_mvgc|>.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] A. B. Barrett, L. Barnett and A. K. Seth, "Multivariate Granger causality
% and generalized variance", _Phys. Rev. E_ 81(4), 2010.
%
% [3] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

% Test data generation

ntrials   = 1;      % number of trials
nobs      = 1000;   % number of observations per trial
fs        = 200;    % sample rate (Hz)
seed      = 0;      % random seed (0 for unseeded)

% VAR model order estimation

moregmode = '';     % model order selection regression mode ('OLS', 'LWR' or empty for default)
mosel     = 'LRT';  % model order selection ('actual', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
momax     = 10;     % maximum model order for model order selection
moalpha   = 0.01;   % significance level for model order selection tests
motstat   = '';     % statistical test for LR model order selection:  'F' for F-test or 'chi2' for chi2 test (default)

% VAR model parameter estimation

regmode   = '';     % VAR model estimation regression mode ('OLS', 'FLS', 'LWR' or empty for default)

% MVGC (time domain) statistical inference

alpha     = 0.05;   % significance level for Granger casuality significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

% MVGC (frequency domain)

fres      = [];     % spectral MVGC frequency resolution (empty for automatic calculation)

%% Generate VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var9_test;
nvars = size(AT,1); % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);
%icfac = 1;
%SQRTV = randn(nvars,icfac*nvars)/sqrt(icfac*nvars); % sample covariance of n x n white noise of length icfac*n
%SIGT  = SQRTV*SQRTV';

% Generate multi-trial VAR time series data with normally distributed residuals
% for specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.

figure(1); clf;
ptic('\n*** tsdata_to_varmo... ');
[moAIC,moBIC,moHQC,moLRT] = tsdata_to_varmo(X,momax,moregmode,moalpha,motstat);
ptoc;

moact = size(AT,3); % actual model order

% Select and report model order.

fprintf('\nModel order selection:');
fprintf('\n\tactual = %d',moact); if strcmpi(mosel,'actual'), morder = moact; fprintf(' *** selected'); end
fprintf('\n\tAIC    = %d',moAIC); if strcmpi(mosel,'AIC'),    morder = moAIC; fprintf(' *** selected'); end
fprintf('\n\tBIC    = %d',moBIC); if strcmpi(mosel,'BIC'),    morder = moBIC; fprintf(' *** selected'); end
fprintf('\n\tHQC    = %d',moHQC); if strcmpi(mosel,'HQC'),    morder = moHQC; fprintf(' *** selected'); end
fprintf('\n\tLRT    = %d',moLRT); if strcmpi(mosel,'LRT'),    morder = moLRT; fprintf(' *** selected'); end
if isnumeric(mosel) && isint(mosel)
    morder = mosel;
    fprintf('\n\tuser   = %d *** selected',morder);
else
    if morder == momax, fprintf(2,' *** WARNING: maximum model order may have been set too low'); end
    if morder == 0,     fprintf(2,' *** ERROR: zero model order! GCs will all be zero!'); return; end
end
fprintf('\n');

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed - bailing out');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Convert VAR to SS

% Convert the VAR model to innovations-form SS model. The state-space dimension
% is morder*nvars. Note that the residuals covariance matrix |SIG| is the same
% for both models, and that - in contrast to the VAR autoregression - the SS
% state-space autoregression is always 1-lag.
%
% _IMPORTANT:_ We also check the SS model for stability, minimum phase and
% symmetric positive-definite residuals covariance matrix. _THIS CHECK SHOULD
% ALWAYS BE PERFORMED!_ - subsequent routines may fail if there are errors here.
% If there are problems with the data (e.g. non-stationarity, colinearity, etc.)
% there's also a good chance they'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong.

fprintf('\n*** var_to_ss\n');
[A1,C,K,ssinfo] = var_to_ss(A,SIG);  % A1 is the "companion matrix" of A. Checks
                                     % are performed, information about the SS
                                     % is returned in 'ssinfo', and (by default)
                                     % a SS diagnostic report is printed out.
% Check for SS errors

assert(~ssinfo.error,'SS error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from SS model.

ptic('*** ss_to_pwcgc... ');
F = ss_to_pwcgc(A1,C,K,SIG);
ptoc;

ptic('*** var_to_pwcgc... ');
FF = var_to_pwcgc(A,SIG);
ptoc;

fprintf(2,'\nmax|F-FF| = %g\n',maxabs(F-FF));

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)])

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% Calculate spectral pairwise-conditional causalities resolution from SS model
% parameters. If not specified, we set the frequency resolution to something
% sensible (based on the spectral radii of the SS model - see ss_info) - we also
% warn if the calculated resolution is very large, as this may cause problems.

if isempty(fres)
    fres = ssinfo.fres;
    if fres > 3000 % adjust to taste
        fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
        istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
    end
end

ptic(sprintf('\n*** ss_to_spwcgc (at frequency resolution = %d)... ',fres));
f = ss_to_spwcgc(A1,C,K,SIG,fres);
ptoc;

ptic(sprintf('\n*** var_to_spwcgc (at frequency resolution = %d)... ',fres));
ff = var_to_spwcgc(A,SIG,fres);
ptoc;

fprintf(2,'\nmax|f-ff| = %g\n',maxabs(f-ff));

% Check for failed spectral GC calculation

assert(~isbad(f,false),'spectral GC calculation failed');

% Plot spectral causal graph.

figure(3); clf;
plot_spw(f,fs);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities. Note that this may occasionally fail if a certain condition of
% Geweke on the VAR parameters is not satisfied: see J. Geweke, J. Am. Stat.
% Assoc. 77(378), 1982 and J. Am. Stat. Assoc. 79(388), 1984.

Fint = bandlimit(f,3); % integrate spectral MVGCs (frequency is dimension 3 of CPSD array (v2.0 - 'smvgc_to_mvgc' deprecated)

fprintf('\n*** GC spectral integral check... ');
res = max(max(abs(F-Fint)./(1+abs(F)+abs(Fint)))); % maximum relative residual
if res < 1e-5
    fprintf('PASS: max relative residual = %.2e\n',res);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',res);
end

%%
% <mvgc_demo.html back to top>
