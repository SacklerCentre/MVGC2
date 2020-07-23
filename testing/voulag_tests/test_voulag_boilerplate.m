%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('T',     'var'), T     = 2;     end % time (secs)
if ~exist('fs',    'var'), fs    = 10000; end % sampling frequency (Hz)
if ~exist('dect',  'var'), dect  = 10;    end % decay time (ms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/fs;
m  = round(T*fs);
T  = m*dt;
t  = (0:m)'*dt;

z = sqrt(dt)*randn(m+1,1);
d = 1-dt/(dect/1000);

x = mvfilter([],d,z')';

gp_qplot(t,x);

[S,f] = periodogram(x,[],[],fs);

gpcmds = [];
gp_qplot(f,log(S),[],gpcmds);
