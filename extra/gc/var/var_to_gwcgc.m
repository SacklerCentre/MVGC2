function F = var_to_gwcgc(A,V,groups)

% Groupwise-conditional GCs
%
% Groups supplied as a cell vector of index vectors.
%
% GC is calculated between each pair of groups (first index is target
% group, second is source group), conditioning on all other variables
% in the system specified by the VAR paremeters A,V.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(groups,n);

F = nan(g);
for b = 1:g
	y = groups{b};
    r = 1:n; r(y) = []; % omit group b
	[~,VR,rep] = vardare(A,V,y,r);
    if sserror(rep,b), continue; end % check DARE report, bail out on error
    for a = 1:g
        if a == b, continue; end
        x = groups{a};
        xr = findin(x,r); % indices of group{a} in r
        F(a,b) = logdet(VR(xr,xr))-logdet(V(x,x));
    end
end
