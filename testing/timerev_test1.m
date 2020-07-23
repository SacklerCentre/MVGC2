%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Minimal VAR parameters

a = 0.8;
b = 0.9;
c = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simple minimal VAR

A = [a c; 0 b];
V = eye(2);

% Calculate autocovariance

G = var_to_autocov(A,V);

q = size(G,3)-1;

fprintf('\nautocovariance lags = %d\n\n',q);

% "Time reversal" - transpose the autocovariance sequence

for k = 1:q+1
	Grev(:,:,k) = G(:,:,k)';
end

% Calculate GCs (Whittle's algorithm) for original and time-reversed autocovs

F    = autocov_to_mvgc(G,   1,2); % original: 2->1
Frev = autocov_to_mvgc(Grev,2,1); % reversed: 1->2

% Check results

fprintf('original: GC(2->1) = %6.4f\n',  F);
fprintf('reversed: GC(1->2) = %6.4f\n\n',Frev);
