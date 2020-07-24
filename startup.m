%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVGC Toolbox "startup" script
%
% Initialise MVGC Toolbox. This file is run automatically if Matlab is started
% in the toolbox root (installation) directory.
%
% You may have to (or want to) customise this script for your computing
% environment.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set toolbox version

global mvgc_version;
mvgc_version.major = 2;
mvgc_version.minor = 0;

fprintf('[mvgc startup] Initialising MVGC toolbox version %d.%d\n', mvgc_version.major, mvgc_version.minor);

% Set paths

% Add mvgc root directory and appropriate subdirectories to path

global mvgc_root;
mvgc_root = fileparts(mfilename('fullpath')); % directory containing this file

% essentials
addpath(mvgc_root);
addpath(fullfile(mvgc_root,'core'));
addpath(genpath(fullfile(mvgc_root,'gc')));
addpath(genpath(fullfile(mvgc_root,'cc')));
addpath(fullfile(mvgc_root,'stats'));
addpath(fullfile(mvgc_root,'utils'));
addpath(fullfile(mvgc_root,'demo'));
addpath(fullfile(mvgc_root,'mex'));
addpath(fullfile(mvgc_root,'docs')); % don't add the 'html' subdirectory!

addpath(fullfile(mvgc_root,'experimental')); % not thouroughly tested - use with care!

% addpath(fullfile(mvgc_root,'GCCA_compat'));    % moved to deprecated
% addpath(fullfile(mvgc_root,'MVGC_v1_compat')); % moved to deprecated
addpath(genpath(fullfile(fullfile(mvgc_root,'deprecated')))); % added v2.0 - make deprecated functions available for the time being

% if ~fexists(@rng) || ~fexists(@randi) % legacy hacks disabled

if true % extra stuff
	addpath(genpath(fullfile(mvgc_root,'extra')));
end

if true % testing (to taste)
	addpath(fullfile(mvgc_root,'testing'));
end

if true % Initialise in-house "Gpmat" Gnuplot/Matlab library if present.
	global gpmat_root;
	gpmat_root = getenv('MATLAB_GPMAT');
	if exist(gpmat_root,'dir') == 7
		cd(gpmat_root);
		startup;
		cd(mvgc_root);
		fprintf('[mvgc startup] Initialised "Gpmat" Matlab Gnuplot API\n');
	else
		fprintf(2,'[mvgc startup] WARNING: couldn''t find "Gpmat" Matlab Gnuplot API\n');
	end
end

if false % Initialise in-house "graphs" GraphViz/Matlab library if present.
	global graphs_root;
	graphs_root = getenv('MATLAB_GRAPHS');
	if exist(graphs_root,'dir') == 7
		cd(graphs_root);
		startup;
		cd(mvgc_root);
		fprintf('[mvgc startup] Initialised "graphs" Matlab GraphViz API\n');
	else
		fprintf(2,'[mvgc startup] WARNING: couldn''t find "graphs" Matlab GraphViz\n');
	end
end

if false % Initialise in-house LZ library
	global LZc_root;
	LZc_root = getenv('MATLAB_LZC');
	if exist(LZc_root,'dir') == 7
		cd(LZc_root);
		startup;
		cd(mvgc_root);
		fprintf('[mvgc startup] Initialised "LZc" Matlab Lempel-Ziv complexity API\n');
	else
		fprintf(2,'[mvgc startup] WARNING: couldn''t find "fLZc" Matlab Lempel-Ziv complexity API\n');
	end
end

% Maintainer

if true % MVGC maintainer
	addpath(fullfile(mvgc_root,'maintainer'));
end

% Check for mex files and set flags appropriately

global have_mvfilter_mex;
have_mvfilter_mex = exist('mvfilter_mex','file') == 3;
if have_mvfilter_mex
	fprintf('[mvgc startup] ''mvfilter_mex'' mex routine available for your platform\n');
else
	fprintf(2,'[mvgc startup] WARNING: no ''mvfilter'' mex file found; please run ''make'' from\n');
	fprintf(2,'[mvgc startup]          the command line in the C subfolder, then ''mextest''\n');
	fprintf(2,'[mvgc startup]          from the Matlab prompt. Meanwhile, a slower scripted\n');
	fprintf(2,'[mvgc startup]          routine will be used.\n');
end

global have_findin_mex;
have_findin_mex = exist('findin_mex','file') == 3;
if have_findin_mex
    fprintf('[mvgc startup] ''findin'' mex routine available for your platform\n');
else
	fprintf(2,'[mvgc startup] WARNING: no ''findin'' mex file found; please run ''make'' from\n');
	fprintf(2,'[mvgc startup]          the command line in the C subfolder, then ''mextest''\n');
	fprintf(2,'[mvgc startup]          from the Matlab prompt. Meanwhile, a slower scripted\n');
	fprintf(2,'[mvgc startup]          routine will be used.\n');
end

% Check for dependencies on some Matlab(R) toolboxes

% Check if we have Statistics toolbox - see if ch2cdf is present v2.0 -
% Statistics Toolbox is more-or-less mandatory (workarounds removed altogether
% as they caused heartache).

if fexists(@chi2cdf)
    fprintf('[mvgc startup] Statistics Toolbox(TM) seems to be present.\n');
else
    fprintf(2,'[mvgc startup] WARNING: Matlab Statistics Toolbox(TM) does not seem to be present.\n');
    fprintf(2,'[mvgc startup]          Some functionality (in particular statistical inference) may not work.\n');
end

% Check if we have Signal Processing toolbox - see if pwelch is present

if fexists(@pwelch)
	fprintf('[mvgc startup] Signal Processing Toolbox(TM) seems to be present.\n');
else
	fprintf(2,'[mvgc startup] WARNING: Matlab Signal Processing Toolbox(TM) does not seem to be present.\n');
	fprintf(2,'[mvgc startup]          Some functionality (in particular spectral estimation) may not work.\n');
end

% Check if we have 'dlyap' from the Control System toolbox v2.0 - only really
% needed in non-essential routines, and the replacement 'lyapslv' works okay, so
% we don't warn, just inform.

if fexists(@dlyap)
	fprintf('[mvgc startup] Control System Toolbox(TM) seems to be present.\n');
else
	addpath(fullfile(mvgc_root,'utils','control'));
	fprintf(2,'[mvgc startup] INFORMATION: Matlab Control System Toolbox(TM) does not seem to be present.\n');
	fprintf(2,'[mvgc startup]              Slightly slower scripted routines will be used (non-essential).\n');
end

% Initialise rng to avoid predictability of sessions

rng_seed(-1); % seed from /dev/urandom (Unix/Mac) else from clock (Windows)

fprintf('[mvgc startup] Random number generator initialised randomly!\n');

% Enable all warnings

warning on all
fprintf('[mvgc startup] All warnings enabled\n');

% Important notes to users

fprintf('[mvgc startup]\n');
fprintf('[mvgc startup] NOTE 1: PLEASE DO NOT ADD THE FULL MVGC HIERARCHY TO YOUR MATLAB SEARCH PATH!\n');
fprintf('[mvgc startup]         Doing so is likely to cause problems. This script has already set up\n');
fprintf('[mvgc startup]         MVGC paths correctly for your Matlab environment.\n');
fprintf('[mvgc startup]\n');
fprintf('[mvgc startup] NOTE 2: It is highly recommended that any single-precision floating-point data\n');
fprintf('[mvgc startup]         be converted to double precision; some routines may be inaccurate or\n');
fprintf('[mvgc startup]         numerically unstable for single-precision input.\n');
fprintf('[mvgc startup]\n');

% Done

fprintf('[mvgc startup] Initialisation complete (you may re-run ''startup'' at any time)\n');
