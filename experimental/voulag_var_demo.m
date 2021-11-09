%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear; fdec = [19; 20; 21]; conx = [2 1 20/1000 1000; 3 2 18/1000 900]; voulag_var_demo
% clear; fdec = [19; 20; 21]; conx = [2 1 10/1000 3000; 3 2 18/1000 2000; 1 2 15/1000 1000]; voulag_var_demo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g     = 0;
vseed = 9071526;
fs    = 10000;
fres  = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(fdec)
c = size(conx,1)

fdec
conx

rstate = rng_seed(vseed);
V = corr_rand(n,g);
rng_restore(rstate);

[A,rho] = voulag_var(fdec,conx,fs);

p = size(A,3)
rho

S = var_to_cpsd(A,eye(n),fres,true);
f = sfreqs(fres,fs);

gp_qplot(f,S,[],'set logs xy');

F = var_to_spwcgc(A,V,fres);
plot_sgc(F,f,[],[],'x11');
