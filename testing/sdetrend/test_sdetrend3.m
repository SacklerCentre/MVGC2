
fs     = 600;    % sample rate (Hz)

fres   = 1000;

%-------------------------------------------------------------------------------

x = X(i,:)';

m = length(x);
fgapp = f + 8*[-1,1]*(fs/m);

t = (0:m)'/fs;

fprintf('calculating sinusoidal fit ... ');
[ff,fgap] = sinufit(x,fs,f,false,1e-12);
fprintf('done\n\n');
disp([f;ff]);

fprintf('calculating MSE ... ');
fE = linspace(fgapp(1),fgapp(2),fres)';
E = sinufit(x,fs,fE,true);
fprintf('done\n');

wl = sprintf('%s\n%s\n%s\n%s',wline(f,2),wline(ff,3),wline(fgap(1),1),wline(fgap(2),1));
gp_qplot(fE,E,[],wl);

function wl = wline(f,ls)

	wl = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 ls %d nohead',f,f,ls);

end
