% Autocovariance sequence for innovations-form SS model. See
%
% L. Barnett and A. K. Seth, Granger casusality for state-space models,
% Phys. Rev. E 91 040101(R), 2015; Supplemental Material, eqs. 1a, 1b.

function [G,q] = ss_to_autocov(A,C,K,V,q)

if nargin < 2 || isempty(q), q = -eps; end % negative q is tolerance

n = ss_parms(A,C,K,V);

if q < 0 % -q is absolute lag
	q = -q;
	alags = false;
else
	alags = true;
end
G = zeros(n,n,q+1);

M = dlyap(A,K*V*K');   % Omega
G(:,:,1) = C*M*C' + V; % Gamma_0
if q == 0, return; end

L = A*M*C' + K*V;      % Lambda_1
G(:,:,2) = C*L;        % Gamma_1

if alags % iterate until ||G_q|| hits floating-point accuracy (relative to G_0)
	tol = eps(maxabs(G(:,:,1)));
	k = 2;
	d = maxabs(G(:,:,k));
	while d > tol && k <= q
		k = k+1;
		L = A*L;           % Lambda_{k-1}
		G(:,:,k) = C*L;    % Gamma_{k-1}
		d = maxabs(G(:,:,k));
	end
	G(:,:,k+1:end) = [];   % truncate uncalculated lags
	q = k-1;
else
	for k = 2:q
		L = A*L;           % Lambda_k
		G(:,:,k+1) = C*L;  % Gamma_k
	end
end
