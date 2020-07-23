%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-space dynamical independence optimisation demo with multiple restarts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default parameters (override on command line)

if ~exist('n',     'var'), n     = 7;     end % microscopic state dimension
if ~exist('r',     'var'), r     = 3*n;   end % hidden state dimension
if ~exist('rhoa',  'var'), rhoa  = 0.9;   end % state-space AR spectral norm (< 1)
if ~exist('g',     'var'), g     = 1.0;   end % residuals multi-information
if ~exist('m',     'var'), m     = 4;     end % macroscopic state dimension
if ~exist('T',     'var'), T     = 10000;     end % macroscopic state dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a random state-space model in innovations form

[A,C,K,rhob] = iss_rand(n,r,rhoa);
V = corr_rand(n,g);

TODO
