T     = 2;       % s

fs    = 10000;  % Hz

mtheta = 2000;   % Hz
wtheta = 0.1;    % width
%gain   = 1;

%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/fs;
m  = round(T*fs);
T  = m/fs;
t  = (0:m)'/fs;

x = sqrt(dt)*randn(m+1,1);
dctx = dct(x);

F = (0:m)'*(fs/m);

mtheta = 2*mtheta;
stheta = wtheta*mtheta;
ga = (mtheta^2)/(stheta^2); % Gamma shape parameter
gb = (stheta^2)/mtheta;     % Gamma scale parameter
p  = gampdf(F,ga,gb);       % Gamma-distributed frequency peak

p = gain*(std(dctx)/max(p))*p;

max(p)
std(dctx)

y = idct(dctx+p);

gp_qplot(t,x);
gp_qplot(F,p);

[S,f] = periodogram(y,[],[],fs);

gpcmds = [];
gp_qplot(f,log(S),[],gpcmds);
