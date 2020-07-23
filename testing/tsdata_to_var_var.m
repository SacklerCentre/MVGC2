function [A,V,E] = tsdata_to_var_var(X,p,dm)

if nargin < 3 || isempty(dm), dm = true;  end

[n,m] = size(X);
assert(p < m,'too many lags');

V = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')
E = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')

if dm, X = demean(X); end % no constant term!

C  = (X*X')/m;

F         = X';
B         = X';
AF        = zeros(n,n*p);
AB        = zeros(n,n*p);
PE        = zeros(n,n*p+n);
PE(:,1:n) = C;
PEF       = C;
PEB       = C;

for k = 1:p

    M  = m-k;
    
    D = (F(k+1:m,:)'*B(1:m-k,:))/M;

    kk = k*n+(1-n:0);
    
    AF(:,kk) = D / PEB;	
    AB(:,kk) = D'/ PEF;

    tmp        = F(k+1:m,:) - B(1:m-k,:)*AF(:,kk).';
    B(1:m-k,:) = B(1:m-k,:) - F(k+1:m,:)*AB(:,kk).';
    F(k+1:m,:) = tmp;

    for r = 1:k-1
        rf = r*n+(1-n:0);
        rb = (k-r)*n+(1-n:0);
        tmp      = AF(:,rf) - AF(:,kk)*AB(:,rb);
        AB(:,rb) = AB(:,rb) - AB(:,kk)*AF(:,rf);
        AF(:,rf) = tmp;
    end;

    PEF = (F(k+1:m,:)'*F(k+1:m,:))/M;
    PEB = (B(1:m-k,:)'*B(1:m-k,:))/M;

    PE(:,k*n+(1:n)) = PEF;        
end;

A = reshape(AF,n,n,p);
