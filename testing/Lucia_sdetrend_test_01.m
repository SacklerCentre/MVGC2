datadir = fullfile(getenv('DATADIR'),'Lucia','sdetrend_test');

datafile = fullfile(datadir,'testdat.mat');

load(datafile);

chans = [31 2 16];
flo = 1;
fhi = 60;

fres = [];
nwin = [];

X = X(chans,:);

leg = num2str(chans','chan %2d');

gpcmds = '\nset logs y\nset grid\nset xlabel "frequency (Hz)"\nset ylabel "power (dB)" rot';

[SX,fX] = tsdata_to_cpsd(X,0,fs,[],[],[],true);

frange = find(fX >= flo & fX <= fhi);
fX = fX(frange);
SX = SX(frange,:);

gp_qplot(fX,SX,leg,['set title "before"' gpcmds]);

ffreq = 3;
fharm = 1:20;
freqs = ffreq*fharm;
wind = [2.0 0.1];

[Y,ts] = sdetrend(X,fs,freqs,wind,ts,1);

[SY,fY] = tsdata_to_cpsd(Y,0,fs,[],[],[],true);

frange = find(fY >= flo & fY <= fhi);
fY = fY(frange);
SY = SY(frange,:);

gp_qplot(fY,SY,leg,['set title "after"' gpcmds]);
