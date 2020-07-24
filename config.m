%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVGC2 Toolbox "startup" script
%
% Configure MVGC2 Toolbox. This script is run by "startup.m".
%
% You may have to (or want to) customise this script for your computing
% environment.
%
% (C) Lionel Barnett and Anil K. Seth, 2020. See file LICENSE in installation
% directory for licensing terms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[mvgc startup] Initialising MVGC2 toolbox\n');

% Set paths

% Add MVGC2 root directory and appropriate subdirectories to path

mvgc2_root = fileparts(mfilename('fullpath')); % directory containing this file

% Optional paths

include_experimental = true;
include_extras       = true;
include_deprecated   = false;
include_testing      = true;
include_maintainer   = false;

gpmat_root = fullfile('~','git','gpmat'); % set to empty to omit
assert(exist(gpmat_root,'dir') == 7,'bad "gpmat" path: ''%s'' does not exist or is not a directory',gpmat_root);

graphs_root = fullfile('~','localrepo','matlab','graphs'); % set to empty to omit
assert(exist(graphs_root,'dir') == 7,'bad "graphs" path: ''%s'' does not exist or is not a directory',graphs_root);

flzc_root = fullfile('~','localrepo','c','fLZc','matlab'); % set to empty to omit
assert(exist(flzc_root,'dir') == 7,'bad "fLZc" path: ''%s'' does not exist or is not a directory',flzc_root);
