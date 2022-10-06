% Autocovariance sequence for innovations-form SS model. See
%
% L. Barnett and A. K. Seth, Granger casusality for state-space models,
% Phys. Rev. E 91 040101(R), 2015; Supplemental Material, eqs. 1a, 1b.

function [G,p] = ss_to_autocov(A,C,K,V,pmax,tol)

% NOTES:
%
% Gamma_k is returned in G(:,:,k+1)
%
% Positive pmax is maximum lags; calculate sequence iteratively until within tolerance.
%
% Negative pmax is absolute number of lags.
%
% if pmax == 0, the covariance matrix Gamma_0 is returned in G.

n = ss_parms(A,C,K,V);

M = dlyap(A,K*V*K'); % Omega
G = C*M*C' + V;      % Gamma_0

if pmax == 0
	p = 0;
	return
end

alags = pmax < 0; % -pmax is absolute number of lags
if alags
	p = -pmax;
	G = cat(3,G,zeros(n,n,p));
end

L = A*M*C' + K*V; % Lambda_1
G(:,:,2) = C*L;   % Gamma_1

if alags
	for k = 2:p
		L = A*L;          % Lambda_k
		G(:,:,k+1) = C*L; % Gamma_k
	end
else
	% iterate until ||Gamma_k|| hits tolerance (default: floating-point accuracy relative to G_0)
	if nargin < 6 || isempty(tol), tol = eps(maxabs(G(:,:,1))); end
	k = 2;
	while maxabs(G(:,:,k)) > tol
		if k > pmax
			fprintf(2,'WARNING: covariance sequence failed to converge (increase max. lags?)\n');
			p = k-1;
			return
		end
		k = k+1;
		L = A*L;        % Lambda_{k-1}
		G(:,:,k) = C*L; % Gamma_{k-1}
	end
	p = k-1;
end
