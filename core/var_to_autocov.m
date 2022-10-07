% Autocovariance sequence for VAR model
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

function [G,q] = var_to_autocov(A,V,qmax,tol)

[n,~,p] = size(A);
pn = p*n;
pn1 = (p-1)*n;

% Associated VAR(1) parameters

A  = reshape(A,n,pn);
A1 = [A; eye(pn1) zeros(pn1,n)]; % "companion matrix" for A
V1 = [V zeros(n,pn1); zeros(pn1,n) zeros(pn1)];

% Solve the Lyapunov equation for the covariance matrix of the associated VAR(1)

G1 = dlyap(A1,V1);

G = reshape(G1(1:n,:),n,n,p);

alags = qmax < 0; % -qmax is absolute number of lags
if alags, q = -qmax; else, q = qmax; end
q1 = q+1;

% We already have p-1 lags; if that's enough, truncate if necessary and return

if q < p
	G = G(:,:,1:q1);
	return
end

% Initialise reverse covariance sequence

R = zeros(pn,n);
for k = 1:p
	R((k-1)*n+1:k*n,:) = G(:,:,p-k+1);
end

% Calculate autocovariances iteratively

if alags % calculate recursively from p lags up to q lags

	G = cat(3,G,zeros(n,n,q1-p)); % pre-allocate
	for k = p+1:q1
		G(:,:,k) = A*R;
		R = [G(:,:,k);R(1:pn1,:)]; % update reverse covariance sequence
	end

else     % calculate recursively from p lags until convergence or maximum lags exceeded

	if nargin < 4 || isempty(tol), tol = eps(maxabs(G(:,:,1))); end
	k = p;
	while maxabs(G(:,:,k)) > tol
		if k > qmax
			fprintf(2,'WARNING: covariance sequence failed to converge (increase max. lags?)\n');
			q = k-1;
			return
		end
		k = k+1;
		G(:,:,k) = A*R;
		R = [G(:,:,k);R(1:pn1,:)]; % update reverse covariance sequence
	end
	q = k-1;

end
