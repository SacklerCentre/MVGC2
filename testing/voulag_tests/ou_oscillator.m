function [x,t,fsmin] = ou_oscillator(fdec,fosc,fs,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fdec      decay frequency (Hz)
% fosc      oscillator frequency (Hz)
% fs        sampling frequency (Hz)
% T         simulation time (secs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/fs;
m  = round(T*fs);
t  = (0:m)'*dt;

nfdec = 2*pi*fdec*dt; % normalised decay frequency (radians)
nfosc = 2*pi*fosc*dt; % normalised oscillator frequency (radians)
fsmin = fs*(nfdec*nfdec+nfosc*nfosc)/(2*nfdec); % minimum sample frequency for stability (Hz)

if fs < fsmin+eps % unstable - bail out
	x = [];
	t = [];
	return
end

z = sqrt(dt)*randn(m+1,2);

x = mvfilter([],[1-nfdec -nfosc; nfosc 1-nfdec],z')';
