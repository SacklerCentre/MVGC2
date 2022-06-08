% Almost an 'idare' drop-in replacement (but info is different)
%
% For efficiency of multiple runs with same size problem, run
%
%	[~,info] = slidare(A,C,Q,R,S);
%	ldwork = info.ldwork;
%
% to get optimal 'ldwork' (workspace size)
%
% Calls 'idare' if slidare_mex not available

function [P,info,V,K,L] = slidare(A,C,Q,R,S,ldwork)

[r, r1] = size(A); assert(r1 == r);
[r1, n] = size(C); assert(r1 == r);
[r1,r2] = size(Q); assert(r1 == r && r2 == r);
if nargin < 4 || isempty(R)
	R = eye(n);
else
	[n1,n2] = size(R); assert(n1 == n && n2 == n);
end
if nargin < 5 || isempty(S)
	S = zeros(r,n);
else
	[r1,n1] = size(S); assert(r1 == r && n1 == n);
end

global have_slidare_mex;
if have_slidare_mex

	if nargin < 6 || isempty(ldwork), ldwork = 0; end
	[P,rep,rcond,ldwork,Lr,Li,Ls] = slidare_mex(A,C,Q,R,S,ldwork);
	if nargout > 1
		info.rep    = rep;
		info.rcond  = rcond;
		info.ldwork = ldwork;
		if nargout > 2
			[W,p] = chol(P);
			if p == 0
				WC = W*C;
				V = WC'*WC+R;
			else
				V = C'*P*C+R;
			end
			if nargout > 3
				K = V\(C'*P*A+S');
				if nargout > 4
					L = [Lr(1:r)+1i*Li(1:r)]./Ls(1:r);
				end
			end
		end
	end

else

	[P,K,L,info1] = idare(A,C,Q,R,S,'noscaling'); % slightly slower, possibly less accurate, but *really* slow with scaling
	if nargout > 1
		info.rep = info1.Report;
		if nargout > 2
			[W,p] = chol(P);
			if p == 0
				WC = W*C;
				V = WC'*WC+R;
			else
				V = C'*P*C+R;
			end
		end
	end

end
