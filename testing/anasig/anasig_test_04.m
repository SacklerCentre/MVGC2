
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fs',    'var'), fs    = 250;       end
if ~exist('tlen',  'var'), tlen  = 15;        end
if ~exist('omega', 'var'), omega = 87;        end
if ~exist('mu',    'var'), mu    = 0;         end
if ~exist('fres',  'var'), fres  = [];        end
if ~exist('fwind', 'var'), fwind = 40;        end
if ~exist('fover', 'var'), fover = [];        end
if ~exist('seed',  'var'), seed  = 0;         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time stamp (seconds)

nobs = round(fs*tlen)+1;
tlen = (nobs-1)/fs;
t = (0:nobs-1)/fs;

% Seed random number generator.

rng_seed(seed);

% Generate sinusoidal time series data

phi = jitter*randn(1,nobs);

X = sin(2*pi*omega*(t-phi)) + mu;

% Calculate analytic signal of bandpass-filtered signal

Z = hilbert(X')';

Z1 = real(Z);
Z2 = imag(Z);

Y = abs(Z);

if true
%	tlim = [5 7];
	mlim = round(tlim*fs);
	mseg = [mlim(1):mlim(2)];
	figure(4);
%	plot(t(mseg)',[X(mseg);Y(mseg)]');
	plot(t(mseg)',[Z1(mseg);Z2(mseg);Y(mseg)]');
	xlim(tlim);
	xlabel('seconds');
%	legend('raw','analytic');
	legend('real','imag','anal');
	return
end

% Remove temporal means

%X = demean(X);
%Y = demean(Y);

[Sx,f] = periodogram(X');
[Sy,f] = periodogram(Y');
f = (fs/(2*pi))*f;

figure(2); clf;
flim = [0 f(end)];
subplot(2,1,1)
plot(f,log(Sx));
xlim(flim);
xlabel('frequency (Hz)');
ylabel('power (dB)');
subplot(2,1,2)
plot(f,log(Sy));
xlim(flim);
xlabel('frequency (Hz)');
ylabel('power (dB)');
