% Autocovariance sequence for innovations-form SS model. See
%
% L. Barnett and A. K. Seth, Granger casusality for state-space models,
% Phys. Rev. E 91 040101(R), 2015; Supplemental Material, eqs. 1a, 1b.
%
% NOTES:
%
% Gamma_k is returned in G(:,:,k+1)
%
% Positive qmax is maximum lags; calculate sequence iteratively until
% within tolerance or maximum lags exceeded.
%
% Negative qmax is absolute number of lags.
%
% if qmax == 0, the covariance matrix Gamma_0 is returned in G.
%
% Default tolerance is machine fp precision relative to Gamma_0.

function [G,q] = ss_to_autocov(A,C,K,V,qmax,tol)

n = ss_parms(A,C,K,V);

M = dlyap(A,K*V*K'); % Omega
G = C*M*C' + V;      % Gamma_0

if qmax == 0
	q = 0;
	return
end

alags = qmax < 0; % -qmax is absolute number of lags
if alags
	q = -qmax;
	G = cat(3,G,zeros(n,n,q)); % pre-allocate
end

L = A*M*C' + K*V; % Lambda_1
G(:,:,2) = C*L;   % Gamma_1

if alags % calculate recursively from 2 lags up to q lags

	for k = 3:q+1
		L = A*L;        % Lambda_k
		G(:,:,k) = C*L; % Gamma_k
	end

else     % calculate recursively from 2 lags until convergence or maximum lags exceeded

	if nargin < 6 || isempty(tol), tol = eps(maxabs(G(:,:,1))); end
	k = 2;
	while maxabs(G(:,:,k)) > tol
		if k > qmax
			fprintf(2,'WARNING: covariance sequence failed to converge (increase max. lags?)\n');
			q = k-1;
			return
		end
		k = k+1;
		L = A*L;        % Lambda_{k-1}
		G(:,:,k) = C*L; % Gamma_{k-1}
	end
	q = k-1;

end
