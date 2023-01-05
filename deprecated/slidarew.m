% Version of slidare with preallocated workspaces (doesn't return eigenvalues)
%
% First run
%
%	[ALPHAR,ALPHAI,BETA,SS,TT,UU,IWORK,DWORK,BWORK] = slidarew_alloc(A,C,Q,R,S);
%
% to get optimal 'ldwork' (workspace size)

function [P,info,V,K] = slidarew(A,C,Q,R,S,ALPHAR,ALPHAI,BETA,SS,TT,UU,IWORK,DWORK,BWORK)

[P,rep,rcond,ldwork] = slidarew_mex(A,C,Q,R,S,ALPHAR,ALPHAI,BETA,SS,TT,UU,IWORK,DWORK,BWORK);

K = (C'*P*C+R)\(C'*P*A+S');

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
		end
	end
end
