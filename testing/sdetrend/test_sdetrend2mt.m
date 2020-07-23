wfac   = 0.4;    % frequency factor (fraction of Nyqvist)
a      = 0.2;    % sinusoid intensity
phi    = 0.7;    % phase
T      = 10;     % data length (seconds)
n      = 3;      % no. of variables
fs     = 600;    % sample rate (Hz)

fres   = 2^14;
winfac = 4;
r      = 0.1;

p   = 5;
rho = 0.8;
w   = 0.5;

%-------------------------------------------------------------------------------

a = [0.1 0.3 0.2]

f = wfac*(fs/2);

m = round(fs*T);
T = m/fs;

fgapp = f + 4*[-1 1]/T;

t = (0:m)'/fs;

s = a.*sin(2*pi*(f*t-phi));

u = var_to_tsdata(var_rand(n,p,rho,w),1,m+1)';

x = u + s;

fprintf('calculating sinusoidal fit ... ');
[ffit,fgap] = sinufitx(x,fs,f,false,1e-12);
fprintf('done\n\n');
disp([f;ffit]);

fprintf('calculating MSE ... ');
fE = linspace(fgapp(1),fgapp(2),fres)';
E = sinufitx(x,fs,fE,true);
fprintf('done\n');

wl = sprintf('%s\n%s\n%s\n%s',wline(f,2),wline(ffit,3),wline(fgap(1),1),wline(fgap(2),1));
gp_qplot(fE,E,[],wl);

y = zeros(size(x));
for i = 1:n
	[y(:,i),~,P(i),Q(i)] = sdetrendx(x(:,i),fs,ffit);
end

yy = sdetrendx(x,fs,ffit);

maxabs(yy-y)

[a hypot(P,Q)]

fprintf('calculating u spectrum ... ');
[Su,freqs,fres] = tsdata_to_cpsd(u',false,fs,winfac,[],fres,true);
fprintf('done\n');

fprintf('calculating x spectrum ... ');
[Sx,freqs,fres] = tsdata_to_cpsd(x',false,fs,winfac,[],fres,true);
fprintf('done\n');

fprintf('calculating y spectrum ... ');
[Sy,freqs,fres] = tsdata_to_cpsd(y',false,fs,winfac,[],fres,true);
fprintf('done\n');

k = find(freqs >= (1-r)*f & freqs <= (1+r)*f);
gp_qplot(freqs(k),[Su(:,k)' Sx(:,k)' Sy(:,k)'],{'original','with sinusoid','sdterend'});

function wl = wline(f,ls)

	wl = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 ls %d nohead',f,f,ls);

end
