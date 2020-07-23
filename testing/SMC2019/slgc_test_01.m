nx   = 3;
ny   = 5;
nz   = 2;
p    = 20;
rho  = 0.9;
w    = 1.0;
g    = 0.5;

K    = [5 8 17];
c    = 2;

m    = 10000;

alpha = 0.05;

gpterm = 'epsl';
wdisp  = '';
mhc    = 'FDR';

stat = 'LR';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx+ny+nz;
u = randperm(n);
x = u(1:nx)
y = u(nx+1:nx+ny)

NK = 1:p; NK(K) = [];
KK = zeros(p,1); KK(K) = 1;

AA = var_rand(n,p,rho,w,wdisp);
AA(x,y,K)  = c*AA(x,y,K);
AA(x,y,NK) = 0;
AA = specnorm(AA,rho);

VV = corr_rand(n,g);

[F2,info] = var_to_slgc(AA,VV,x,y,p,verb);
info

U = varfima_to_tsdata(AA,[],[],VV,m);

[~,stats] = tsdata_to_slgc(U,x,y,p,verb);
F1  = stats.(stat).tstat;
F1c = stats.(stat).icdf(1-alpha);
F1p = 1-stats.(stat).cdf(F1);
F1s = significance(F1p,alpha,mhc);

[F1 F2]

gpstem = fullfile(tempdir,mfilename);
gp_write(gpstem,[(1:p)' F1 F2]);
gp = gp_open(gpstem,gpterm,[1 2]);
bgap = 0.2;
ymfac = 0.05;
ygfac = 0.01;
ymax = (1+ymfac)*max([F1;F2]);
fprintf(gp,'datfile = "%s.dat"\n',gpstem);
fprintf(gp,'set xlabel "lag"\n');
fprintf(gp,'set ylabel "GC" norot\n');
fprintf(gp,'set xtics 1\n');
fprintf(gp,'set xr [%g:%g]\n',(1-bgap)/2,p+(1+bgap)/2);
fprintf(gp,'set boxwidth %g relative\n',1-bgap);
fprintf(gp,'set style fill transparent solid 0.2\n');
fprintf(gp,'unset key\n');
fprintf(gp,'set grid y\n');
fprintf(gp,'set yr [0:%g]\n',ymax);

fprintf(gp,'\nset multiplot layout 1,2\n');

fprintf(gp,'set title "F (sample)"\n');
fprintf(gp,'unset arrow; set arrow from graph 0,first %g to graph 1,first %g front nohead lw 2 lc rgb "red"\n',F1c,F1c);
fprintf(gp,'plot datfile u 1:2 w boxes\n');

fprintf(gp,'set title "F (analytic)"\n');
fprintf(gp,'plot datfile u 1:3 w boxes\n');

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,2);
