%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVGC2 Toolbox "config" script
%
% Configure MVGC2 Toolbox. This script is run by "startup.m".
%
% This is the default configuration script; to retain local configuration options,
% copy this file to "local_config.m" and customise.
%
% (C) Lionel Barnett and Anil K. Seth, 2020. See file LICENSE in installation
% directory for licensing terms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configure paths, etc.

% Add MVGC2 root directory and appropriate subdirectories to path

mvgc2_root = fileparts(mfilename('fullpath')); % directory containing this file

% Optional paths (defaults - override on command line)

if ~exist('include_experimental', 'var'),include_experimental = true;  end
if ~exist('include_extras',       'var'),include_extras       = true;  end
if ~exist('include_deprecated',   'var'),include_deprecated   = false; end
if ~exist('include_testing',      'var'),include_testing      = true;  end
if ~exist('include_maintainer',   'var'),include_maintainer   = false; end

if ~exist('gpmat_root', 'var')
	gpmat_root = '';  % set empty to omit
end

if ~exist('graphs_root', 'var')
	graphs_root = ''; % set empty to omit
end

if ~exist('flzc_root', 'var')
	flzc_root = '';   % set empty to omit
end
