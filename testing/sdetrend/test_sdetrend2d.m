f  = 50;
fs = 600;


fres   = 2^14;
winfac = 2;
r      = 0.3;

%-------------------------------------------------------------------------------
%
% E.g.,
%
% trial = 1; chan = 20; x = (10^12)*double(dat{trial}(chan,:)');
%-------------------------------------------------------------------------------

m = length(x);
T = m/fs;

fgapp = f + 4*[-1 1]/T;

t = (0:m)'/fs;

fprintf('calculating sinusoidal fit ... ');
[ffit,fgap] = sinufit(x,fs,f,false,1e-12);
fprintf('done\n\n');
disp([f;ffit]);

fprintf('calculating MSE ... ');
fE = linspace(fgapp(1),fgapp(2),fres)';
E = sinufit(x,fs,fE,true);
fprintf('done\n');

wl = sprintf('%s\n%s\n%s\n%s',wline(f,2),wline(ffit,3),wline(fgap(1),1),wline(fgap(2),1));
gp_qplot(fE,E,[],wl);

[y,~,P,Q] = sdetrendx(x,fs,ffit);
hypot(P,Q)

fprintf('calculating x spectrum ... ');
[Sx,freqs,fres] = tsdata_to_cpsd(x',false,fs,winfac,[],fres);
Sx = squeeze(Sx);
fprintf('done\n');

fprintf('calculating y spectrum ... ');
[Sy,freqs,fres] = tsdata_to_cpsd(y',false,fs,winfac,[],fres);
Sy = squeeze(Sy);
fprintf('done\n');

k = find(freqs >= (1-r)*f & freqs <= (1+r)*f);
gp_qplot(freqs(k),[Sx(k) Sy(k)],{'original','sdterend'});

function wl = wline(f,ls)

	wl = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 ls %d nohead',f,f,ls);

end
