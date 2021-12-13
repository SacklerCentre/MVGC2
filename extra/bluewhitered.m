function bwr = bluewhitered(sparms,sfun,len,fignum)

% Make blue-white-red colormap with reasonable defaults.
% 'sfun' specifies a "stretch function", taking parameters
% 'sparms'; it must either be a string specifying a built-in
% stretch function, or a function handle.

if nargin < 2 || isempty(sfun), sfun = 'pow'; end
if nargin < 3 || isempty(len),  len  = 128;   end

if ischar(sfun)
	nsfun = sfun;
	switch lower(sfun)
	case 'pow'
		if nargin < 1 || isempty(sparms), sparms = 0.5; end
		sfun = @(x) realpow(x,sparms);
	case 'cospow'
		if nargin < 1 || isempty(sparms), sparms = 0.4; end
		sfun = @(x) 1-cos((pi/2)*realpow(x,sparms));
	case 'exp'
		if nargin < 1 || isempty(sparms), sparms = 0.4; end
		a = 1/(1-exp(-1/(2*sparms^2)));
		sfun = @(x) a*(1-exp(-((x/sparms).^2)/2));
	otherwise
		error('Unknown "stretch function"');
	end
else
	nsfun = 'user-supplied';
	assert(isa(sfun,'function_handle'),'Must supply a "stretch function" as a function handle');
end

% Have the stretch function

mask = sfun((0:len-1)'/(len-1));
red = [ones(len,1) mask mask];
blu = [mask mask ones(len,1)];
bwr = [blu;flipud(red)];

% For testing

if nargin > 3 || ~isempty(fignum)
	figure(fignum);
	clf;
	x = linspace(0,1,500)';
	y = sfun(x);
	plot([-flipud(y);y],[-flipud(x);x]);
	grid on;
	xlabel('scale');
	ylabel('intensity');
	title(sprintf('Stretch function : %s\n',nsfun));
	colormap(bwr);
	colorbar;
	caxis([-1,1]);
end
