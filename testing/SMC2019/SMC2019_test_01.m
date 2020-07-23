%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supply: m, LRstats

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

rstate = rng_seed(useed);
U = varfima_to_tsdata(A,[],[],V,m);
rng_restore(rstate);

[~,stats] = tsdata_to_pwslgc(U,p);

[FF,info] = var_to_pwslgc(A,V,p);
info

F  = stats.(stat).tstat;
Fp = 1-stats.(stat).cdf(F);
Fc = stats.(stat).icdf(1-alpha/nmh);
Fs = significance(Fp,alpha,'Bonferroni');

conest  = [];
for i = 1:n
	for j = 1:n
		if j == i, continue; end
		for u = 1:p
			if Fs(u,i,j), conest = [conest;[i,j,u]]; end
		end
	end
end

diag(A(:,:,1))'
con
conest

k = 0;
data = (1:p)';
for i = 1:n
	for j = 1:n
		if j == i, continue; end
		k = k+1;
		data = [data F(:,i,j)];
	end
end

gpstem = fullfile(tempdir,mfilename);
gp_write(gpstem,data);
gp = gp_open(gpstem,gpterm,[1.8 0.6]);
bgap = 0.1;
ymfac = 0.1;
ygfac = 0.01;
fgap = 0.2;
xlp = 0.25;
ylp = 1.02;
tol = sqrt(eps);
fprintf(gp,'datfile = "%s.dat"\n',gpstem);
fprintf(gp,'set xlabel "lag"\n');
%fprintf(gp,'set ylabel "GC" rot\n');
fprintf(gp,'unset ylabel\n');
fprintf(gp,'set xtics 5 scale 0\n');
fprintf(gp,'set xr [0:%g]\n',p+1);
fprintf(gp,'set boxwidth %g relative\n',1-bgap);
fprintf(gp,'set style fill transparent solid 0.2\n');
fprintf(gp,'unset ytics\n');
fprintf(gp,'unset key\n');
fprintf(gp,'set linestyle 9 lt 1 lc "black" lw 6\n');
fprintf(gp,'set grid y\n');
fprintf(gp,'set multiplot title "%s"layout %d,%d\n',ptitle,n,n);

samescale = 1;

if samescale
	maxF = max(FF(:));
	maxFF = max(F(:));
	ymax = (1+ymfac)*max(maxF,maxFF);
	fprintf(gp,'set yr [0:%g]\n',ymax);
end

k = 0;
for i = 1:n
	for j = 1:n
		if j == i
			fprintf(gp,'set multiplot next\n');
			continue
		end
		k = k+1;
		fprintf(gp,'set title "%s"\n',sprintf('$%d \\\\to %d$',j,i));
		if ~samescale
			ymax = (1+ymfac)*max([max(data(:,k+1)) Fc max(max(FF(k,:,:)))]);
			fprintf(gp,'set yr [0:%g]\n',ymax);
		end
		fprintf(gp,'set arrow from first 0,first %g to graph 1,first %g front nohead lw 2 lc rgb "red"\n',Fc,Fc);
		for v = 1:p
			if FF(v,i,j) > tol
				fprintf(gp,'set arrow from first %g,first %g to first %g,first %g front nohead ls 9\n',v-0.5-fgap,FF(v,i,j),v+0.5+fgap,FF(v,i,j));
			end
		end

		fprintf(gp,'plot datfile u 1:%d w boxes not\n',1+k);
		fprintf(gp,'unset arrow\n');
	end
end

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,2);
