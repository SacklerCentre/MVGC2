% All times in milliseconds!

n   = 5;
T   = 5*1000;
%dt = 0.01;

a = [4 6 5 3 5];

c = [ ...
	1 2  8 30; ...
	1 4  7 25; ...
	1 5  9 22; ...
	2 4 10 40; ...
	3 2  9 25; ...
    3 4  5 18; ...
	5 1  8 20 ...
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = ceil(T/dt)+1;
T = dt*(m-1)/1000; % seconds
t = (0:m)/1000;    % seconds

x = randn(n,m+1);

fprintf('\nscripted version : '); tic
ym = sdde(a,c,x,dt,false);
et = toc; fprintf('%8.6f seconds\n',et);

fprintf('mex version      : '); tic
yx = sdde(a,c,x,dt,true);
et = toc; fprintf('%8.6f seconds\n',et);

fprintf('\nmax. absolute difference = %g\n\n',maxabs(ym-yx));

gp_qplot(t',yx',[],sprintf('set xr [0:%g]',T),[],[1 Inf])
