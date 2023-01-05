%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVGC2 Toolbox default configuration script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the default MVGC configuration script. To set user configuration options,
% copy this file to "config.m" in this directory and edit to taste.
%
% The configuration script is run by 'startup'; any of these configuration options
% may be overriden on the command line before 'startup' is called.
%
% (C) Lionel Barnett and Anil K. Seth, 2020. See file LICENSE in installation
% directory for licensing terms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal options

if ~exist('include_experimental', 'var'), include_experimental = true;  end
if ~exist('include_extras',       'var'), include_extras       = true;  end
if ~exist('include_deprecated',   'var'), include_deprecated   = false; end
if ~exist('include_testing',      'var'), include_testing      = false; end

% Optional external library paths - set to empty to omit from initialisation
%
% Examples
%
% if ~exist('gpmat_path',  'var'), gpmat_path = getenv('GPMAT_PATH'); end
% if ~exist('gvmat_path',  'var'), gvmat_path = getenv('GVMAT_PATH'); end

if ~exist('gpmat_path', 'var'), gpmat_path = ''; end
if ~exist('gvmat_path', 'var'), gvmat_path = ''; end
