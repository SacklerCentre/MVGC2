function [KR,VR,rep] = vardare(VARA,V,y,r)

% Solve DARE for SS model derived from reduced VAR (see var2iss)

ny   = length(y);
nr   = length(r);
p    = size(VARA,3);
pny  = p*ny;
pny1 = pny-ny;

A = [reshape(VARA(y,y,:),ny,pny); eye(pny1) zeros(pny1,ny)];
C = reshape(VARA(r,y,:),nr,pny);
Q = [V(y,y) zeros(ny,pny1); zeros(pny1,pny)];
S = [V(y,r); zeros(pny1,nr)];
R = V(r,r);
[KR,VR,rep] = mdare(A,C,Q,R,S);
