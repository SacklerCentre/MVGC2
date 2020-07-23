
fres = 2048;
fs = 200;

SAR = cpsd2autospec(varfima_to_cpsd([],[],[d r true ],1,fres));
SMA = cpsd2autospec(varfima_to_cpsd([],[],[d r false],1,fres));

f = sfreqs(fres,fs);

STH = realpow(2-2*cos(f/fs),-d);
%STH = realpow(f/fs,-2*d);

gp_qplot(f,[SAR' SMA' STH],{'AR','MA','theoretical'},'set logs xy')
