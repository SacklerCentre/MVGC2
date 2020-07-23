function [A,V,E] = tsdata_to_var_LWR_ref(X,p,dm)

if nargin < 3 || isempty(dm), dm = true;  end

[n,m] = size(X);
assert(p < m,'too many lags');

V = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')
E = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')

if dm, X = demean(X); end % no constant term!

X = X';

PE  = covm(X);

F = X;
B = X;
PEF = PE;
PEB = PE;
for k = 1:p,
        D	= covm(F(k+1:m,:),B(1:m-k,:));

        AF(:,k*n+(1-n:0)) = D / PEB;	
        AB(:,k*n+(1-n:0)) = D'/ PEF;

        tmp        = F(k+1:m,:) - B(1:m-k,:)*AF(:,k*n+(1-n:0)).';
        B(1:m-k,:) = B(1:m-k,:) - F(k+1:m,:)*AB(:,k*n+(1-n:0)).';
        F(k+1:m,:) = tmp;

        for r = 1:k-1,
                tmp = AF(:,r*n+(1-n:0))   - AF(:,k*n+(1-n:0))*AB(:,(k-r)*n+(1-n:0));
                AB(:,(k-r)*n+(1-n:0)) = AB(:,(k-r)*n+(1-n:0)) - AB(:,k*n+(1-n:0))*AF(:,r*n+(1-n:0));
                AF(:,r*n+(1-n:0))   = tmp;
        end;

        PEF = covm(F(k+1:m,:),F(k+1:m,:));

        PEB = covm(B(1:m-k,:),B(1:m-k,:));

        PE(:,k*n+(1:n)) = PEF;        
end;

A = reshape(AF,n,n,p);

function C = covm(x,y)

if nargin < 2
    C = x'*x;
else
    C = x'*y;
end
C = C/(size(x,1));
