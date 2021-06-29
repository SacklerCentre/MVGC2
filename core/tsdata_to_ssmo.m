function [mosvc,rmax] = tsdata_to_ssmo(X,pf,plotm)

% Estimate SS model porder using Bauer's Singular Value Criterion (SVC; D. Bauer, Automatic 37, 2001)
%
% X       - observation process time series
% pf      - past/future horizons for canonical correlations
% plotm   - empty: don't plot, integer: Matlab plot to figure n (if zero, use next); string: Gnuplot terminal (may be empty)
%
% mosvc   - Bauer's Singular Value Criterion (SVC) optimal model order
% rma     - maximum possible model order given supplied past/future horizons
%
% The past/future horizons pf may be supplied as a 2-vector [p,f] or a scalar p
% = f = pf. Bauer (D. Bauer, Automatic 37, 2001) recommends setting p = f = 2*p,
% where p is the optimal VAR model order for the observation process X according
% to Aikaike's Information Criterion (AIC).

[n,m,N] = size(X);

assert(all(isint(pf(:))),'past/future horizon must be a 2-vector or a scalar positive integer');
if isscalar(pf)
    p = pf;    f = pf;
elseif isvector(pf) && length(pf) == 2
    p = pf(1); f = pf(2);
else
    error('past/future horizon must be a 2-vector or a scalar positive integer');
end
assert(p+f < m,'past/future horizon too large (or not enough data)');
rmax = n*min(p,f);

X = demean(X); % no constant term (don't normalise!)

mp  = m-p;
mp1 = mp+1;
mf  = m-f;
mh  = mp1-f; % m-p-f+1

M  = N*mp;
M1 = N*mp1;
Mh = N*mh;

Xf = zeros(n,f,mh,N);
for k = 1:f
    Xf(:,k,:,:) = X(:,p+k:mf+k,:);
end
Xf = reshape(Xf,n*f,Mh);

XP = zeros(n,p,mp1,N);
for k = 0:p-1
    XP(:,k+1,:,:) = X(:,p-k:m-k,:);
end
Xp = reshape(XP(:,:,1:mh,:),n*p,Mh);
XP = reshape(XP,n*p,M1);

[Wf,cholp] = chol((Xf*Xf')/Mh,'lower');
assert(cholp == 0,'forward weight matrix not positive definite');

[Wp,cholp] = chol((Xp*Xp')/Mh,'lower');
assert(cholp == 0,'backward weight matrix not positive definite');

BETA = Xf/Xp; % 'OH' estimate: regress future on past
assert(all(isfinite(BETA(:))),'subspace regression failed');

[~,S] = svd(Wf\BETA*Wp); % SVD of CCA-weighted OH estimate

sval = diag(S);       % the singular values
df   = 2*n*(1:rmax)'; % number of free parameters (Hannan & Deistler, see also Bauer 2001) ... or rmax*rmax+2*n*rmax ???
svc  = -log(1-[sval(2:rmax);0]) + df*(log(Mh)/Mh); % Bauer's Singular Value Criterion

morder = (0:rmax)';
[~,idx] = min(svc); mosvc = morder(idx);

if ~isempty(plotm) % we're going to plot

	mo = (1:rmax)';

	if mosvc == rmax, wsvc = '*'; else wsvc = ''; end

	gap = 0.05;
	ssvc = gap+(1-gap)*(svc-min(svc))/(max(svc)-min(svc));

	if ischar(plotm) % Gnuplot

		gpname = 'sssvc';
		gpstem = fullfile(tempdir,gpname);
		gpdat = [mo ssvc sval];
		gp_write(gpstem,gpdat);

		gp = gp_open(gpstem,plotm,[Inf,0.6]);

		fprintf(gp,'datfile = "%s.dat"\n',gpname);

		fprintf(gp,'\nset grid\n');
		fprintf(gp,'set xr[0:%g]\n',rmax);
		fprintf(gp,'set xlabel "SS dimension"\n');

		fprintf(gp,'\nset multiplot title "SS SVC model order selection (CCA, max = %d)\\\n" layout 2,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n',rmax);

		fprintf(gp,'\nset title "Singular value criterion (SVC)"\n');
		fprintf(gp,'set ytics 0.2\n');
		fprintf(gp,'set yr[0:1.05]\n');
		fprintf(gp,'set ylabel "SVC (scaled)"\n');
		fprintf(gp,'set key top right Left rev\n');
		fprintf(gp,'plot \\\n');
		fprintf(gp,'datfile u 1:2 w linespoints pt 6 ps 1.4 t "SVC (opt = %2d%c)"\n',mosvc,wsvc);

		fprintf(gp,'\nset title "Singular values"\n');
		fprintf(gp,'unset logs y\n');
		fprintf(gp,'set ytics auto format ''%% h''\n');
		fprintf(gp,'set yr [0:*]\n');
		fprintf(gp,'set ytics 0.2 nomirror\n');
		fprintf(gp,'set ylabel "singular value"\n');
		fprintf(gp,'plot datfile u 1:3 w boxes fs solid 0.25 not\n');

		fprintf(gp,'\nunset multiplot\n');

		gp_close(gp,gpstem,plotm);

	else % Matlab

		if plotm == 0, figure; else, figure(plotm); end; clf;

		xlims = [0 rmax];

		subplot(2,1,1);
		plot(mo,ssvc,'o-');
		grid on
		title('Singular value criterion (SVC)');
		ylabel('SVC (scaled)');
		xlabel('SS dimension');
		legend(sprintf('SVC (opt = %d%c)',mosvc,wsvc));
		xlim(xlims);
		ylim([0 1+gap]);

		subplot(2,1,2); % SVC
		bar(mo,sval,1.01,'FaceColor',[0.65 0.75 1]);
		xlim(xlims);
		xlabel('SS dimension');
		ylabel('singular value');
		title('singular values');

		axes('Units','Normal');
		h = title(sprintf('SS SVC model order selection (CCA, max = %d)\n\n',rmax),'FontSize',13);
		set(gca,'visible','off')
		set(h,'visible','on')
	end
end
