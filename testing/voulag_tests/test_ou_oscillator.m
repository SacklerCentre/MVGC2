%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fdec',  'var'), fdec  = 10;    end % decay frequency (Hz)
if ~exist('fosc',  'var'), fosc =  100;   end % oscillator frequency (Hz)
if ~exist('fs',    'var'), fs    = 10000; end % sampling frequency (Hz)
if ~exist('T',     'var'), T     = 2;     end % time (secs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,t,fsmin] = ou_oscillator(fdec,fosc,fs,T);

assert(fs > fsmin,'Unstable: set sample frequency > %g Hz',fsmin);

x = x(:,1);

gp_qplot(t,x);

[S,f] = periodogram(x,[],[],fs);

gpcmds = sprintf('set arrow from first %g,graph 0 to first %g,graph 1 nohead ls 8 lw 2',fosc,fosc);
gpcmds = [gpcmds '\nset xr [1:500]'];
% gpcmds = [gpcmds '\nset logs x'];
gp_qplot(f,log(S),[],gpcmds);
