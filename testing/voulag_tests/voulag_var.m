function [A,rho] = voulag_var(fdec,conx,fs,rho)

n = length(fdec);
assert(isvector(fdec),'Node decays must be a vector');

fdecn = 2*pi*fdec/fs; % normalised node decay frequencies (radians)

[c,cols] = size(conx);
assert(ismatrix(conx) && cols == 4,'Connectivity table must be a 4-column matrix');

itarget = conx(:,1);
isource = conx(:,2);
ilag    = round(fs*conx(:,3))+1; % causal feedback lag
fcause  = 2*pi*conx(:,4)/fs;     % normalised causal feedback frequencies (radians)

assert(all(itarget >= 1 & itarget <= n),'Some target indices out of range');
assert(all(isource >= 1 & isource <= n),'Some source indices out of range');

A = zeros(n,n,max(ilag));

fdec
fdecn
1-fdecn

for i = 1:n
	A(i,i,1) = 1-fdecn(i);
end

for k = 1:c
	A(itarget(k),isource(k),ilag(k)) = fcause(k);
end

if nargin < 4 || isempty(rho)
	if nargout > 1
		rho = specnorm(A);
	end
else
	A = specnorm(A,rho);
end
