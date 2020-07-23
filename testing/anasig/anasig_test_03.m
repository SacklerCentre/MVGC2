
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fs',    'var'), fs    = 250;       end
if ~exist('tlen',  'var'), tlen  = 100;       end
if ~exist('nvars', 'var'), nvars = 5;         end
if ~exist('ssmo',  'var'), ssmo  = 11;        end
if ~exist('mu',    'var'), mu    = 4.0;       end
if ~exist('rhoa',  'var'), rhoa  = 0.9;       end
if ~exist('rmi',   'var'), rmi   = 0.5;       end
if ~exist('ford',  'var'), ford  = 20;        end
if ~exist('frip',  'var'), frip  = 1;         end
if ~exist('fattn',  'var'),fattn = 40;        end
if ~exist('fband', 'var'), fband = [80 100];  end
if ~exist('fres',  'var'), fres  = [];        end
if ~exist('fwind', 'var'), fwind = 100;       end
if ~exist('fover', 'var'), fover = [];        end
if ~exist('seed',  'var'), seed  = 0;         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time stamp (seconds)

nobs = round(fs*tlen)+1;
tlen = (nobs-1)/fs;
t = (0:nobs-1)'/fs;

% Seed random number generator.

rng_seed(seed);

% Generate random SS parameters

[A,C,K] = iss_rand(nvars,ssmo,rhoa);
V = corr_rand(nvars,rmi);
ssinfo = ss_info(A,C,K,V);
assert(~ssinfo.error,'SS error(s) found - bailing out');

% Generate SS time series data

Xraw = ss_to_tsdata(A,C,K,V,nobs) + mu;

% Band-pass raw signal
nfband = fband/(fs/2);

%[z,p,k] = butter(ford,nfband); sos = zp2sos(z,p,k);
[z,p,k] = ellip(ford,frip,fattn,nfband); sos = zp2sos(z,p,k);

[isstable(sos) isminphase(sos)]

%figure(1); fvtool(sos);
figure(1); freqz(sos); ylim([-100 10]);

for i =1:nvars
	Xbpf(i,:) = filtfilt(sos,1,Xraw(i,:)')';
end

% Calculate analytic signal of bandpass-filtered signal

%Xana = abs(hilbert(Xbpf')');
Xana = envelope(Xbpf')';

if false
%	i = 3;
%	tlim = [5 7];
	mlim = round(tlim*fs);
	mseg = [mlim(1):mlim(2)];
	figure(4);
	subplot(2,1,1);
	plot(t(mseg),Xraw(i,mseg)');
	xlim(tlim);
	xlabel('seconds');
	legend('raw');
	subplot(2,1,2);
	plot(t(mseg),[Xbpf(i,mseg); Xana(i,mseg)]');
	xlim(tlim);
	legend({'filtered','envelope'});
	xlabel('seconds');
end

% Remove temporal means
%%{
Xraw = demean(Xraw);
Xbpf = demean(Xbpf);
Xana = demean(Xana);
%%}

% Frequency resolution

if isempty(fres)
	fres = 2^nextpow2(ssinfo.acdec);
%	fres = 2^nextpow2(nobs);
	fprintf('Using frequency resolution %d\n',fres);
end
if fres > 10000 % adjust to taste
	fprintf(2,'WARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

f = sfreqs(fres,fs);

Sraw = cpsd2apsd(ss_to_cpsd(A,C,K,V,fres));

Srae = tsdata_to_cpsd(Xraw,fs,fwind,fover,fres,true);

Sbpf = tsdata_to_cpsd(Xbpf,fs,fwind,fover,fres,true);

Sana = tsdata_to_cpsd(Xana,fs,fwind,fover,fres,true);

figure(2);
flim = [0 f(end)];
for i = 1:nvars
	subplot(nvars,2,2*i-1);
	plot(f,log([Sraw(:,i) Srae(:,i)]));
	xlim(flim);
	xlabel('frequency (Hz)');
	ylabel('power (dB)');
	legend({'actual','estimated'},'location','northeastoutside');
	subplot(nvars,2,2*i);
	plot(f,log([Sbpf(:,i) Sana(:,i)]));
	xlim(flim);
	xlabel('frequency (Hz)');
	ylabel('power (dB)');
	legend({'filtered','envelope'},'location','northeastoutside');
	line([fband(1) fband(1)],ylim,'HandleVisibility','off');
	line([fband(2) fband(2)],ylim,'HandleVisibility','off');
end
