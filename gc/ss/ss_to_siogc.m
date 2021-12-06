function f = ss_to_siogc(A,C,K,V,inout)

% in/out spectral GCs per variable

[n,r] = ss_parms(A,C,K,V);

if strcmpi(inout,'in')
	gcin = true;
elseif strcmpi(inout,'out')
	gcin = false;
else
	error('in/out parameter must be ''in'' or ''out''');
end

TODO
