%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate some iEEG-like data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All times in seconds, frequencies in Hz!

if ~exist('d',    'var'), d    = 100;    end % duration (seconds)
if ~exist('tx',   'var'), tx   = 10;     end % relaxation time (seconds)
if ~exist('fsim', 'var'), fsim = 100000; end % simulation frequency (Hz)
if ~exist('fs',   'var'), fs   = 1000;   end % "recording" frequency (Hz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = [4 6 5 3 5];

c = [ ...
	1 2  8 30; ...
	1 3  7 25; ...
	1 5  9 22; ...
	2 4 10 40; ...
	3 2  9 25; ...
    3 4  5 18; ...
	5 1  8 20 ...
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(a);

d    = d+tx;
msim = ceil(d*fsim)+1;
d    = (msim-1)/fsim;
tsim = (0:msim-1)/fsim;

% Simulate

fprintf('simulating SDDE ... '); tic
x = sdde(a,c,randn(n,msim),1/fsim,true);
et = toc; fprintf('%8.6f seconds\n',et);

% Truncate

mx   = ceil(tx*fsim);
x    = x(:,mx+1:end);

% Downsample

ds = round(fsim/fs);
fs = fsim/ds;

X = downsample(x,ds);
m = size(X,2);
t = (0:m-1)/fs;

% Plot

if false
	msim = size(x,2);
	tsim = (0:msim-1)/fsim;
	gp_qplot(tsim',x',[],sprintf('set xr [0:%g]',tsim(end)),[],[1 Inf])
end

gp_qplot(t',X',[],sprintf('set xr [0:%g]',t(end)),[],[1 Inf])
