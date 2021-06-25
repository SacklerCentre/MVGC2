%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SS model generation parameters

nvars     = 5;      % number of variables
ssmo      = 9;      % SS model order
rhoa      = 0.95;   % AR spectral radius
rmi       = 0.5;    % residuals log-generalised correlation (multi-information)
                    % g = -log|R|. g = 0 yields zero correlation,g = [] is uniform random
                    % on space of correlation matrices
seed      = 0;      % random smodel eed (0 for unseeded)

% MVGC (frequency domain)

fres      = [];     % spectral MVGC frequency resolution (empty for automatic calculation)

% Spectral factorisarion (Wilson's algorithm) parameters

sftol     = [];     % numerical tolerance (empty for default 1e-7)
sfmaxi    = [];     % maximum iterations  (empty for default 100)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed random number generator.

rng_seed(seed);

% Generate random SS parameters in innovations form
%
% For real data, you would replace this by estimating
% an SS model from the data

[A,C,K] = iss_rand(nvars,ssmo,rhoa);
V = corr_rand(nvars,rmi);

% Report information on the generated SS model and check for errors.

ssinfo = ss_info(A,C,K,V);
assert(~ssinfo.error,'SS error(s) found - bailing out');

% Calculate a reasonable frequency resolution (if not specified)

if isempty(fres)
	maxfres = 2^14;
    fres = calc_fres(rhoa);
	if fres > maxfres
		fprintf(2,'\nWARNING: esitmated frequency resolution %d exceeds maximum; setting to %d\n' ,fres,maxfres);
		fres = maxfres;
	else
		fprintf('\nUsing frequency resolution %d\n',fres);
	end
end

% Calculate CPSD from SS parameters

S = ss_to_cpsd(A,C,K,V,fres);

% Granger causality calculation: time domain, SS

ptic('*** ss_to_pwcgc...');
F1 = ss_to_pwcgc(A,C,K,V);
ptoc;
assert(~isbad(F1,false),'SS/GC estimation failed');

% Granger causality calculation: time domain, CPSD

ptic('\n*** cpsd_to_pwcgc... \n');
F2 = cpsd_to_pwcgc(S,sftol,sfmaxi);
ptoc;
assert(~isbad(F2,false),'CPSD/GC estimation failed');

% Time domain relative error

RE = maxabs(F2(:)-F1(:))/maxabs(F1(:));
fprintf('\nTime domain: SS vs. CPSD relative error = %e\n',RE);

% Granger causality calculation: frequency domain, SS

ptic('\n*** ss_to_spwcgc...');
f1 = ss_to_spwcgc(A,C,K,V,fres);
ptoc;
assert(~isbad(f1,false),'spectral SS/GC calculation failed');

ptic('\n*** cpsd_to_spwcgc...\n');
f2 = cpsd_to_spwcgc(S,sftol,sfmaxi);
ptoc;
assert(~isbad(f2,false),'spectral CPSD/GC calculation failed');

% Time domain relative error

re = maxabs(f2(:)-f1(:))/maxabs(f1(:));
fprintf('\nFrequency domain: SS vs. CPSD relative error = %e\n',RE);

% Check that spectral causalities average (integrate) to time-domain
% causalities. Note that this may occasionally fail if a certain condition
% on the VAR parameters is not satisfied (Geweke 1982).

Fint1 = bandlimit(f1,3);

fprintf('\n*** GC spectral integral check (SS)... ');
rr1 = abs(F1-Fint1)./(1+abs(F1)+abs(Fint1)); % relative residuals
mrr1 = max(rr1(:));                       % maximum relative residual
if mrr1 < 1e-5
    fprintf('PASS: max relative residual = %.2e\n',mrr1);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',mrr1);
end

Fint2 = bandlimit(f2,3);

fprintf('\n*** GC spectral integral check (CPSD)... ');
rr2 = abs(F2-Fint2)./(1+abs(F2)+abs(Fint2)); % relative residuals
mrr2 = max(rr2(:));                       % maximum relative residual
if mrr2 < 1e-5
    fprintf('PASS: max relative residual = %.2e\n',mrr2);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',mrr2);
end
