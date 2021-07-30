%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVGC2 Toolbox start-up script
%
% Initialise MVGC2 Toolbox. This script is run automatically if Matlab is started
% in the toolbox root (installation) directory.
%
% Set configuration options in "config.m".
%
% (C) Lionel Barnett and Anil K. Seth, 2020. See file LICENSE in installation
% directory for licensing terms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mvgc2_root
mvgc2_root = fileparts(mfilename('fullpath')); % directory containing this file

fprintf('[MVGC2 startup] Initialising MVGC2 toolbox\n');

% Set configuration options: look first for config.m in this directory, else run default

% if exist('config.m','file') == 2 is unreliable, as there might already be
% a config.m on the PATH, which 'exist' will pick up!!! Thus we simply look
% up the actual filename in root directory listing.

if isempty(dir('config.m'))
	fprintf('[MVGC2 startup] Setting default configuration options\n');
	config_default;
else
	fprintf('[MVGC2 startup] Setting user configuration options\n');
	config;
end

% Set paths

% Add MVGC2 root directory and appropriate subdirectories to path

addpath(mvgc2_root);
addpath(fullfile(mvgc2_root,'core'));
addpath(genpath(fullfile(mvgc2_root,'gc')));
addpath(genpath(fullfile(mvgc2_root,'cc')));
addpath(fullfile(mvgc2_root,'stats'));
addpath(fullfile(mvgc2_root,'utils'));
addpath(fullfile(mvgc2_root,'demo'));
addpath(fullfile(mvgc2_root,'mex'));
addpath(fullfile(mvgc2_root,'docs')); % don't add the 'html' subdirectory!

if include_experimental
	addpath(fullfile(mvgc2_root,'experimental')); % not thouroughly tested - use with care!
end
clear include_experimental

if include_deprecated
	addpath(genpath(fullfile(fullfile(mvgc2_root,'deprecated'))));
end
clear include_deprecated

if include_extras
	addpath(genpath(fullfile(mvgc2_root,'extra')));
end
clear include_extras

if include_testing
	addpath(genpath(fullfile(mvgc2_root,'testing')));
end
clear include_testing

if include_maintainer % MVGC maintainer
	addpath(fullfile(mvgc2_root,'maintainer'));
end
clear include_maintainer

fprintf('[MVGC2 startup] Internal paths set\n');

if ~isempty(gpmat_path) % Initialise in-house "gpmat" Gnuplot/Matlab library if present.
	assert(exist(gpmat_path,'dir') == 7,'bad "gpmat" path: ''%s'' does not exist or is not a directory',gpmat_path);
	cd(gpmat_path);
	startup;
	cd(mvgc2_root);
	fprintf('[MVGC2 startup] Initialised "gpmat" Matlab Gnuplot API\n');
end
clear gpmat_path

if ~isempty(gvmat_path) % Initialise in-house "gvmat" GraphViz/Matlab library if present.
	assert(exist(gvmat_path,'dir') == 7,'bad "gvmat" path: ''%s'' does not exist or is not a directory',gvmat_path);
	cd(gvmat_path);
	startup;
	cd(mvgc2_root);
	fprintf('[MVGC2 startup] Initialised "gvmat" Matlab GraphViz API\n');
end
clear gvmat_path

if ~isempty(flzc_path) % Initialise in-house LZ library
	assert(exist(flzc_path,'dir') == 7,'bad "fLZc" path: ''%s'' does not exist or is not a directory',flzc_path);
	cd(flzc_path);
	startup;
	cd(mvgc2_root);
	fprintf('[MVGC2 startup] Initialised "fLZc" Matlab Lempel-Ziv complexity API\n');
end
clear flzc_path

% Check for mex files and set flags appropriately

global have_mvfilter_mex;
have_mvfilter_mex = exist('mvfilter_mex','file') == 3;
if have_mvfilter_mex
	fprintf('[MVGC2 startup] ''mvfilter_mex'' mex routine available for your platform\n');
else
	fprintf(2,'[MVGC2 startup] WARNING: no ''mvfilter'' mex file found; please run ''make'' from\n');
	fprintf(2,'[MVGC2 startup]          the command line in the C subfolder, then ''mextest''\n');
	fprintf(2,'[MVGC2 startup]          from the Matlab prompt. Meanwhile, a slower scripted\n');
	fprintf(2,'[MVGC2 startup]          routine will be used.\n');
end

global have_findin_mex;
have_findin_mex = exist('findin_mex','file') == 3;
if have_findin_mex
    fprintf('[MVGC2 startup] ''findin'' mex routine available for your platform\n');
else
	fprintf(2,'[MVGC2 startup] WARNING: no ''findin'' mex file found; please run ''make'' from\n');
	fprintf(2,'[MVGC2 startup]          the command line in the C subfolder, then ''mextest''\n');
	fprintf(2,'[MVGC2 startup]          from the Matlab prompt. Meanwhile, a slower scripted\n');
	fprintf(2,'[MVGC2 startup]          routine will be used.\n');
end

% Check for dependencies on some Matlab(R) toolboxes

% Check if we have Statistics toolbox - see if ch2cdf is present v2.0 -
% Statistics Toolbox is more-or-less mandatory (workarounds removed altogether
% as they caused heartache).

if fexists(@chi2cdf)
    fprintf('[MVGC2 startup] Statistics Toolbox(TM) seems to be present.\n');
else
    fprintf(2,'[MVGC2 startup] WARNING: Matlab Statistics Toolbox(TM) does not seem to be present.\n');
    fprintf(2,'[MVGC2 startup]          Some functionality (in particular statistical inference) may not work.\n');
end

% Check if we have Signal Processing toolbox - see if pwelch is present

if fexists(@pwelch)
	fprintf('[MVGC2 startup] Signal Processing Toolbox(TM) seems to be present.\n');
else
	fprintf(2,'[MVGC2 startup] WARNING: Matlab Signal Processing Toolbox(TM) does not seem to be present.\n');
	fprintf(2,'[MVGC2 startup]          Some functionality (in particular spectral estimation) may not work.\n');
end

% Check if we have 'dlyap' from the Control System toolbox v2.0 - only really
% needed in non-essential routines, and the replacement 'lyapslv' works okay, so
% we don't warn, just inform.

if fexists(@dlyap)
	fprintf('[MVGC2 startup] Control System Toolbox(TM) seems to be present.\n');
else
	addpath(fullfile(mvgc2_root,'utils','control'));
	fprintf(2,'[MVGC2 startup] INFORMATION: Matlab Control System Toolbox(TM) does not seem to be present.\n');
	fprintf(2,'[MVGC2 startup]              Slightly slower scripted routines will be used (non-essential).\n');
end

% Initialise rng to avoid predictability of sessions

rng_seed(-1); % seed from /dev/urandom (Unix/Mac) else from clock (Windows)

fprintf('[MVGC2 startup] Random number generator initialised randomly!\n');

% Enable all warnings

warning on all
fprintf('[MVGC2 startup] All warnings enabled\n');

% Important notes to users

fprintf('[MVGC2 startup]\n');
fprintf('[MVGC2 startup] NOTE 1: PLEASE DO NOT ADD THE FULL MVGC HIERARCHY TO YOUR MATLAB SEARCH PATH!\n');
fprintf('[MVGC2 startup]         Doing so is likely to cause problems. This script has already set up\n');
fprintf('[MVGC2 startup]         MVGC2 paths correctly for your Matlab environment.\n');
fprintf('[MVGC2 startup]\n');
fprintf('[MVGC2 startup] NOTE 2: It is highly recommended that any single-precision floating-point data\n');
fprintf('[MVGC2 startup]         be converted to double precision; some routines may be inaccurate or\n');
fprintf('[MVGC2 startup]         numerically unstable for single-precision input.\n');
fprintf('[MVGC2 startup]\n');

% Done

fprintf('[MVGC2 startup] Initialisation complete (you may re-run ''startup'' at any time)\n');
