function lam = vardec(A,plotm,pres)

% VAR coefficients sequence decay: returns the fitted exponential
% decay factor for the given VAR coefficients sequence.
%
% A       - VAR coefficients sequence (3D array, last index is lag)
% plotm   - empty: don't plot, integer: Matlab plot to figure n (if zero, use next); string: Gnuplot terminal (may be empty)
% pres    - resolution for plotting exponential fit (default: 1000)
%
% lam     - VAR coefficients exponential decay factor (empty for "natural" decay)

p = size(A,3);
a = zeros(p,1);
for k = 1:p
	a(k) = norm(A(:,:,k));
end

parms = [(1:p)' ones(p,1)]\log(a); % fit to exponential - linear regression (OLS)
lam = -parms(1);                   % fitted exponential decay factor

if nargin > 1 && ~isempty(plotm) % we're going to plot
	if nargin < 3 || isempty(pres), pres = 1000; end
	mu =  parms(2); % fitted exponential offset
	pp = (1:p)';
	pf = linspace(0,p+1,pres)';
	af = exp(mu-lam*pf);
	if ischar(plotm) % plotm is gpterm
		gpname = 'vardec';
		gpstem = fullfile(tempdir,gpname);
		gpdat = {[pp a],[pf af]};
		gp_write(gpstem,gpdat);
		gp = gp_open(gpstem,plotm);
		fprintf(gp,'datfile = "%s.dat"\n',gpname);
		fprintf(gp,'set title "\\\\large VAR coefficients magnitude vs lag"\n');
		fprintf(gp,'set key box top right Left rev spacing 1.3 height 0.5\n');
		fprintf(gp,'set xr [0:%d]\n',p+1);
		fprintf(gp,'set xlabel "lags ($k$)"\n');
		fprintf(gp,'set ylabel "VAR coefficients magnitude" rot\n');
		fprintf(gp,'set key top right Left rev\n');
		fprintf(gp,'set grid\n');
		fprintf(gp,'plot datfile i 0 u 1:2 w points pt 7 ps 0.8 t "$\\\\Vert A_k \\\\Vert$", ');
		fprintf(gp,'datfile i 1 u 1:2 w lines lw 3 t "exponential fit ($\\\\lambda = %g$)\n',lam);
		gp_close(gp,gpstem,plotm);
	elseif plotm
		if plotm == 0, figure; else, figure(plotm); end; clf;
		plot(pp,a,'.','markersize',12);
		hold on
		plot(pf,af);
		hold off
		xlabel('lags ($k$)','Interpreter','latex');
		y=ylabel('VAR coefficients magnitude','Interpreter','latex');
		set(y, 'position', get(y,'position')-[0.2,0,0]);
		legend('$\Vert A_k\Vert$',sprintf('exponential fit (%c\\lambda = %g%c)','$',lam,'$'),'Interpreter','latex');
		xlim([0,p+1]);
		title(sprintf('VAR coefficients magnitude vs lag\n'),'fontsize',11,'fontweight','normal');
		grid on
	end
end
