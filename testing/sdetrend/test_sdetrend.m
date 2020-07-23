w      = 100;     % frequency (Hz)
a      = 0.1;    % sinusoid intensity
phi    = 0.3;    % phase
T      = 200;     % data length (seconds)
fs     = 600;    % sample rate (Hz)

fres   = 2048;

%-------------------------------------------------------------------------------

e = 4/fs;

fgap = w+[-e e];
fgapp = w + 8*[-e e];

m = round(fs*T);
T = m/fs;

t = (0:m)'/fs;

s = a*sin(2*pi*(w*t-phi));

x = randn(size(t)) + s;

fprintf('calculating spectrum ... ');
[S,f,fres] = tsdata_to_cpsd(x',1,fs);
S = squeeze(S);
fprintf('done\n');

fprintf('calculating sinusoidal fit ... ');
ww = sinufit(x,fs,fgap,false,1e-12);
fprintf('done\n\n');
disp([w;ww]);

fprintf('calculating MSE ... ');
ff = linspace(fgapp(1),fgapp(2),fres)';
E = sinufit(x,fs,ff,true);
fprintf('done\n');

fprintf('detrending ... ');
xx = sdetrendx(x,fs,w);
fprintf('done\n');

fprintf('calculating spectrum ... ');
[SS,f,fres] = tsdata_to_cpsd(xx',1,fs);
SS = squeeze(SS);
fprintf('done\n\n');

k = find(f >= 90 & f <= 110);
%gp_qplot(f(k),[SS(k) S(k)]);

wl = sprintf('%s\n%s\n%s\n%s',wline(w,2),wline(ww,3),wline(fgap(1),1),wline(fgap(2),1));
gp_qplot(ff,E,[],wl);


function wl = wline(w,ls)

	wl = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 ls %d nohead',w,w,ls);

end
