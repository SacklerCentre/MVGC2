% Almost an 'idare' drop-in replacement (but info is different)
%
% For efficiency of multiple runs with same size problem, run
%
%	[~,~,~,info] = slidare(A,C,Q,R,S);
%	ldwork = info.ldwork;
%
% to get optimal 'ldwork' (workspace size)

function [P,K,L,info] = slidare(A,C,Q,R,S,ldwork)

if nargin < 6 || isempty(ldwork), ldwork = 0; end

[P,rep,rcond,ldwork,Lr,Li,Ls] = slidare_mex(A,C,Q,R,S,ldwork);

K = (C'*P*C+R)\(C'*P*A+S');
m = size(P,1);
L = [Lr(1:m)+1i*Li(1:m)]./Ls(1:m);

info.rep    = rep;
info.rcond  = rcond;
info.ldwork = ldwork;
