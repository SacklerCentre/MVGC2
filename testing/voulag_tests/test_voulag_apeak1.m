T     = 2;       % s

fs    = 10000;  % Hz

mtheta = 2000;   % Hz
%wtheta = 0.1;    % width
%gain   = 1;

%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/fs;
m  = round(T*fs);
T  = m/fs;
t  = (0:m)'/fs;

y = pwnoise(m+1,fs,mtheta,wtheta,gain);

gp_qplot(t,y);
%gp_qplot(F,p);

[S,f] = periodogram(y,[],[],fs);

gpcmds = [];
gp_qplot(f,log(S),[],gpcmds);
