%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVGC2 Toolbox user configuration script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal options

if ~exist('include_experimental', 'var'), include_experimental = true;  end
if ~exist('include_extras',       'var'), include_extras       = true;  end
if ~exist('include_deprecated',   'var'), include_deprecated   = false; end
if ~exist('include_testing',      'var'), include_testing      = true;  end
if ~exist('include_maintainer',   'var'), include_maintainer   = true;  end

% Optional external library paths - set to empty to omit from initialisation

if ~exist('gpmat_path', 'var'), gpmat_path = getenv('GPMAT_PATH'); end
if ~exist('gvmat_path', 'var'), gvmat_path = getenv('GVMAT_PATH'); end
if ~exist('flzc_path',  'var'), flzc_path  = getenv('FLZC_PATH');  end
