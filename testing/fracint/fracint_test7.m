%-------------------------------------------------------------------------------

p         = 6;      % model order
rho       = 0.8;    % spectral radius
wvar      = 0.5;    % var coefficients decay weighting factor

fir       = 1000;

m         = 10000;
mtrunc    = 10*m;

S         = 100;

d         = 0.4;
q         = 1000;

%-------------------------------------------------------------------------------

%{
A = var_rand(1,p,rho,wvar);
h = zeros(S,1);
for s = 1:S
	fprintf('sample %3d of %3d\n',s,S);
	%X = varfima_to_tsdata(A,[],{d fir true},1,m,1,mtrunc);
	X = gendiff(randn(1,m+mtrunc),[d fir false]);
	X = X(mtrunc+1:mtrunc+m);
%gp_qplot((1:m)',X)
%return
%	h(s) = Hurst_est(X',q);
%	X = randn(m,1);
	h(s) = genhurst(X,1);
end

[d mean(h) std(h)]
%}

X = gendiff(randn(1,m+mtrunc),[d fir false]);
X = X(mtrunc+1:mtrunc+m);
h = genhurst(X,2);
fmin = 0.1; fmax = 10;
X = gendiff(randn(1,m+mtrunc),[d fir true]);
[S,f,fres] = tsdata_to_cpsd(X,true,200,[],[],-4,true,1); % autspectra only
fidx = find(f>=fmin&f<=fmax);
[a,b] = fitline(log(f(fidx)'),log(S(:,fidx)));
line1 = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead',fmin,fmin);
line2 = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead',fmax,fmax);
gpcmds = sprintf('set logs xy\n%s\n%s',line1,line2);
ff = [nan(fidx(1)-1,1);f(fidx);nan(fres-fidx(end)+1,1)];
Sf = exp(a*log(ff)'+b);
gp_qplot(f,[S' Sf'],{sprintf('d = %6.4f',d),sprintf('d = %6.4f, H =  %6.4f',-a/2,h)},gpcmds,'png',Inf,16);

function [a,b] = fitline(x,y)

	n = size(x,2);
	xm  = sum(x,2)/n;
	ym  = sum(y,2)/n;
	xxm = sum(x.*x,2)/n;
	xym = sum(x.*y,2)/n;
	a   = (xym-xm.*ym)./(xxm-xm.*xm);
	b   = ym-xm.*a;

end
