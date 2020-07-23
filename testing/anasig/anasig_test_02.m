
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fs',    'var'), fs    = 250;       end
if ~exist('tlen',  'var'), tlen  = 10;        end
if ~exist('varmo', 'var'), varmo = 6;         end
if ~exist('mu',    'var'), mu    = 4.0;       end
if ~exist('rho',   'var'), rho   = 0.9;       end
if ~exist('wvar',  'var'), wvar  = 0.5;       end
if ~exist('rmi',   'var'), rmi   = 0.5;       end
if ~exist('mosel', 'var'), mosel = 'HQC';     end
if ~exist('momax', 'var'), momax = 3*varmo;   end
if ~exist('fres',  'var'), fres  = [];        end
if ~exist('ford',  'var'), ford  = 4;         end
if ~exist('fband', 'var'), fband = [80 120];  end
if ~exist('seed',  'var'), seed  = 0;         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time stamp (seconds)

nobs = round(fs*tlen)+1;
tlen = (nobs-1)/fs;
t = (0:nobs-1)'/fs;

nobs

% Seed random number generator.

rng_seed(seed);

% Generate random VAR coefficients for test network.

AA = var_rand(tnet5,varmo,rho,wvar);
nvars = size(AA,1); % number of variables

% Generate random residuals covariance (in fact correlation) matrix.

VV = corr_rand(nvars,rmi);
infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

% Generate VAR time series data

X = varfima_to_tsdata(AA,[],[],VV,nobs) + mu;

% Band-pass raw signal

[b,a] = butter(ford,fband/(fs/2)); % bandpass
%fvtool(b,a);
for i =1:nvars
	XX(i,:) = filtfilt(b,a,X(i,:)')';
end

% Calculate analytic signal of bandpass-filtered signal

Y = abs(hilbert(XX')');

if true
%	i = 3;
%	tlim = [5 7];
	mlim = round(tlim*fs);
	mseg = [mlim(1):mlim(2)];
	figure(4);
	subplot(2,1,1);
	plot(t(mseg),X(i,mseg)');
	xlim(tlim);
	xlabel('seconds')
	legend('raw')
	subplot(2,1,2);
	plot(t(mseg),[XX(i,mseg); Y(i,mseg)]');
	xlim(tlim);
	legend({'filtered','envelope'})
	xlabel('seconds')
	return
end

% Remove temporal means

X = demean(X);
Y = demean(Y);

% Model order estimation

[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(Y,momax,'LWR',[],[],1);

% Select and report VAR model order.

morder = moselect(sprintf('VAR model order selection (max = %d)',momax), mosel,'ACT',varmo,'AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
assert(morder > 0,'selected zero model order! GCs will all be zero!');
if morder >= momax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end

% Estimate VAR model of selected order from data.

[A,V] = tsdata_to_var(Y,morder,'LWR');
assert(~isbad(A),'VAR estimation failed - bailing out');
info = var_info(A,V);
assert(~info.error,'VAR error(s) found - bailing out');

% Granger causality calculation: time domain on bandpassed analytic signal

F = var_to_pwcgc(A,V,Y,'LWR');
assert(~isbad(F,false),'GC estimation failed');

% Granger causality (actual) calculation: frequency domain, then bandlimit

if isempty(fres)
    fres = 2^nextpow2(max(info.acdec,infoo.acdec)); % alternatively, fres = 2^nextpow2(nobs);
	fprintf('Using frequency resolution %d\n',fres);
end
if fres > 10000 % adjust to taste
	fprintf(2,'WARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

ff = var_to_spwcgc(AA,VV,fres);
assert(~isbad(ff,false),'spectral GC estimation failed');
FF = bandlimit(ff,3,fs,fband); % integrate spectral MVGCs across frequency band

% Plot GCs

maxF = 1.1*max(nanmax(F(:),nanmax(FF(:))));
plot_gc({FF,F},{'PWCGC (actual)','PWCGC (analytic)'},[],[maxF maxF],2);

% Relative error

relerr = (FF-F)/maxabs(FF);
fprintf('\nRelative error =\n');
disp(relerr);
