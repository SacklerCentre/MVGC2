%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supply: LRstats, kmax

if ~exist('alpha','var'), alpha = 0.05; end
if ~exist('useed','var'), useed = 0;    end

gpterm = 'epsl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SMC2019_testvar5;

if LRstats
	stat = 'LR';
	ptitle = 'Likelihood-ratio statistics\n\n';
else
	stat = 'F';
	ptitle = 'F-statistics\n\n';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VARA = A;
[A,C,K] = var_to_ss(VARA,V,true);

DFMS = nan(kmax,ncons);
AFMS = nan(kmax,ncons);
for k = 1:ncons
	fprintf('MS: target %d of %d\n',k,ncons);
	i = con(k,1);
	j = con(k,2);
	DFMS(:,k) = ss_to_msgc(A,C,K,V,i,j,kmax);
	AFMS(:,k) = cumsum(DFMS(:,k));
end

DFFF = nan(kmax,ncons);
AFFF = nan(kmax,ncons);
for k = 1:ncons
	fprintf('FF: target %d of %d\n',k,ncons);
	i = con(k,1);
	j = con(k,2);
	AFFF(:,k) = ss_to_ffgc(A,C,K,V,i,j,kmax);
	DFFF(:,k) = diff([0;AFFF(:,k)]);
end

diag(VARA(:,:,1))'
con

data = (1:kmax)';
for k = 1:ncons
	data = [data DFMS(:,k) AFMS(:,k) DFFF(:,k) AFFF(:,k)];
end

gpstem = fullfile(tempdir,mfilename);
gp_write(gpstem,data);
gp = gp_open(gpstem,gpterm,[1.8 0.6]);
ymfac = 0.03;
yt = [0.1 0.1 0.2 0.1 0.2]'
fprintf(gp,'datfile = "%s.dat"\n',gpstem);
fprintf(gp,'set xlabel "lag"\n');
%fprintf(gp,'set ylabel "GC" rot\n');
fprintf(gp,'unset ylabel\n');
fprintf(gp,'set xtics %g scale 0\n',2);
fprintf(gp,'set xr [1:%g]\n',32);
fprintf(gp,'unset key\n');
fprintf(gp,'set linestyle 1 lt 1 lc rgb "#0000AA" lw 2\n');
fprintf(gp,'set linestyle 2 lt 1 lc rgb "#AA0000" lw 2\n');
fprintf(gp,'set linestyle 3 lt 1 lc rgb "#009900" lw 2\n');
fprintf(gp,'set linestyle 4 lt 1 lc rgb "#8A2BE2" lw 2\n');
fprintf(gp,'set linestyle 9 lt 1 lc rgb "#080808" lw 1\n');
fprintf(gp,'set grid\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set multiplot title "%s"layout %d,%d\n',ptitle,3,2);

samescale = 0;

if samescale
	%ymax = (1+ymfac)*nanmax([DFMS(:);AFMS(:);DFFF(:);AFFF(:)]);
	ymax = (1+ymfac)*nanmax([DFMS(:);AFFF(:)]);
	fprintf(gp,'set yr [0:%g]\n',ymax);
end

for k = 1:ncons
	i = con(k,1);
	j = con(k,2);
	%fprintf(gp,'\nset title "%s"\n',sprintf('$%d \\\\to %d$',j,i));
	fprintf(gp,'\nset label "%s" at graph 0.05,graph 0.95\n',sprintf('$%d \\\\to %d$',j,i));
	fprintf(gp,'set ytics %g\n',yt(k));
	if k == ncons
		fprintf(gp,'set key at graph 1.2,graph 1 left Left rev\n');
	end
	if ~samescale
		ymax = (1+ymfac)*nanmax([DFMS(:,k);AFFF(:,k)]);
		fprintf(gp,'set yr [0:%g]\n',ymax);
	end
	fprintf(gp,'set arrow from first %g,graph 0 to first %g,graph 1 nohead ls 9\n',con(k,3),con(k,3));
	fprintf(gp,'plot \\\n');
	fprintf(gp,'datfile u 1:%d w lines ls 1 t "multi-step GC",\\\n',2+4*(k-1));
	%fprintf(gp,'datfile u 1:%d w lines ls 3 not,\\\n',3+4*(k-1));
	%fprintf(gp,'datfile u 1:%d w lines ls 4 not,\\\n',4+4*(k-1));
	fprintf(gp,'datfile u 1:%d w lines ls 2 t "full-future GC"\n',5+4*(k-1));
	fprintf(gp,'unset label\n');
	fprintf(gp,'unset arrow\n');
end

fprintf(gp,'\nunset multiplot\n');
gp_close(gp,gpstem,gpterm,2);
