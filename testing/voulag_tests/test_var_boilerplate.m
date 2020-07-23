m  = 10000;
fs = 200;  % Hz
%	a  = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%

t  = (0:m)'/fs;


z = randn(1,m+1);

x = mvfilter([],a,z)';
%x = var_to_tsdata(a,1,m+1)';

gp_qplot(t,x);

%[S,f,nwobs,noobs,nwins] = tsdata_to_cpsd(x',fs,window,[],fres,true);
[S,f] = periodogram(x,[],[],fs);

gpcmds = [];
gp_qplot(f,S,[],gpcmds);

% gpcmds = 'set logs x';
% gp_qplot(f,log(S),[],gpcmds);
