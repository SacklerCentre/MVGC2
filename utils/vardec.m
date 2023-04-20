function [A,a] = vardec(A,lam,plotm)

% VAR coefficients sequence decay
%
% [lam,a] = vardec(A) returns the fitted exponential decay factor
%           for the given VAR coefficients sequence A.
%
% [A,a]   = vardec(A,lam) returns the coefficient sequence A weighted
%           so that the decay factor = lam.
%
% Optionally, the norms at each lag are returned in a.
%
% Note that the second calling form will NOT in general preserve the
% spectral radius of the coefficients sequence.

p = size(A,3);
a = zeros(p,1);
for k = 1:p
	a(k) = norm(A(:,:,k));
end

noweight = nargin < 2 || isempty(lam);

if noweight % no decay factor supplied; return fitted decay factor

	parms = [(1:p)' ones(p,1)]\log(a); % fit to exponential (linear regression = OLS)
	A = -parms(1);                     % this is the fitted decay factor lam

else        % weight the AR coefficients so decay = lam

	aa = exp(-lam*(1:p)');
	for k = 1:p
		A(:,:,k) = (aa(k)/a(k))*A(:,:,k);
	end

end

if nargin > 2 && ~isempty(plotm) % we're going to plot
	if ischar(plotm) % plotm is gpterm
		if noweight
			gp_qplot((1:p)',a,[],'set xlabel "k"\nset ylabel "|A|"\nset xtics 1\nunset key\nset grid\nset title "VAR coefficients: magnitude vs lag"',plotm);
		else
			gp_qplot((1:p)',[a aa],{'old','new'},'set xlabel "k"\nset ylabel "|A|"\nset xtics 1\nset key top right\nset grid\nset title "VAR coefficients: magnitude vs lag"',plotm);
		end
	elseif plotm
		if plotm == 0, figure; else, figure(plotm); end; clf;
		if noweight
			plot((1:p)',a);
		else
			plot((1:p)',[a aa]);
			legend('old','new');
		end
		xlabel('$k$','Interpreter','latex');
		y=ylabel('$\Vert A_k \Vert$','rot',0,'Interpreter','latex');
		set(y, 'position', get(y,'position')-[0.2,0,0]);
		title('VAR coefficients: magnitude vs lag');
		grid on
	end
end
