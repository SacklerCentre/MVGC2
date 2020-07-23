%{
clear; tau = [2 2.2]; gam = [-0.5 0.8]; power_law_spec_test_01
%}

T  = 200;
fs = 500;

dt = 1/fs;
m  = floor(fs*T)+1;
T  = (m-1)*dt;
A  = [1-dt/tau(1) dt/gam(1); dt/gam(2) 1-dt/tau(2)];
V  = dt*eye(2);
t  = linspace(0,T,m)';

rho = var_specrad(A);

fprintf('\nsamples = %d\n',  m  );
fprintf(  'rho     = %g\n\n',rho);

% Generate an (approximate) Ornstein-Uhlenbeck time series

x = varfima_to_tsdata(A,[],[],V,m)';

gpcmds = sprintf('set title "OU process: fs = %g Hz, tau = (%g,%g), gamma = (%g,%g) sec"\nset xzeroaxis lt -1\nset grid',fs,tau(1),tau(2),gam(1),gam(2));
gp_qplot(t,x,[],gpcmds,[],[1 Inf]);


[S,f,fres] = tsdata_to_cpsd(x',false,fs,[],[],[],true,true);
SS = var_to_cpsd(A,V,fres);
ST(:,1) = SS(1,1,:);
ST(:,2) = SS(2,2,:);
gpcmds = sprintf('set title "PSD: fs = %g Hz, tau = (%g,%g), gamma = (%g,%g) sec"\nset logs xy\nset grid',fs,tau(1),tau(2),gam(1),gam(2));
%gp_qplot(f,[S ST],[],gpcmds);
gp_mplot({[f S(:,1) ST(:,1)],[f S(:,2) ST(:,1)]},[],[],gpcmds,[],[1 Inf]);
