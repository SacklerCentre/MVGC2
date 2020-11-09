% Autocovariance sequence for innovations-form SS model. See
%
% L. Barnett and A. K. Seth, Granger casusality for state-space models,
% Phys. Rev. E 91 040101(R), 2015; Supplemental Material, eqs. 1a, 1b.

function [G,q] = ss_to_autocov(A,C,K,V,q)

n = ss_parms(A,C,K,V);

M = dlyap(A,K*V*K');   % Omega
G = zeros(n,n,q+1);
G(:,:,1) = C*M*C' + V; % Gamma_0
L = A*M*C' + K*V;      % Lambda_1
G(:,:,2) = C*L;        % Gamma_1
for k = 2:q
	L = A*L;           % Lambda_{k+1}
	G(:,:,k+1) = C*L;  % Gamma_{k+1}
end
