% Allocate workspace for slidarew

function [ALPHAR,ALPHAI,BETA,SS,TT,UU,IWORK,DWORK,BWORK] = slidarew_alloc(A,C,Q,R,S)


[~,info] = slidare(A,C,Q,R,S); % get optimal 'ldwork' for problem

n      = size(A,1);
m      = size(R,1);
n2     = 2*n;
lds    = n2+m;
ldt    = n2+m;
ldu    = n2;
liwork = max(m,n2);
ldwork = info.ldwork;
lbwork = n2;

ALPHAR = zeros(n2,1);
ALPHAI = zeros(n2,1);
BETA   = zeros(n2,1);
SS     = zeros(lds,lds);
TT     = zeros(ldt,n2);
UU     = zeros(ldu,n2);
IWORK  = int64(zeros(liwork,1));
DWORK  = zeros(ldwork,1);
BWORK  = int64(zeros(ldwork,1));
