function bwr = bluewhitered(sparms,sfun,len)

% Make blue-white-red colormap with reasonable defaults.
% 'sfun' specifies a "stretch function", taking parameters
% 'sparms'; it must either be a string specifying a built-in
% stretch function, or a function handle.

if nargin < 2 || isempty(sfun), sfun = 'pow'; end
if nargin < 3 || isempty(len),  len  = 128;   end

if ischar(sfun)
	switch lower(sfun)
	case 'pow'
		if nargin < 1 || isempty(sparms), sparms = 0.75; end
		sfun = @(x) realpow(x,sparms);
	case 'cospow'
		if nargin < 1 || isempty(sparms), sparms = 0.40; end
		sfun = @(x) 1-cos((pi/2)*realpow(x,sparms));
	case 'exp'
		if nargin < 1 || isempty(sparms), sparms = 0.50; end
		a = 1/(1-exp(-1/(2*sparms^2)));
		sfun = @(x) a*(1-exp(-((x/sparms).^2)/2));
	otherwise
		error('Unknown "stretch function"');
	end
else
	assert(isa(sfun,'function_handle'),'Must supply a "stretch function" as a function handle');
end
mask = sfun((0:len-1)'/(len-1));
red = [ones(len,1) mask mask];
blu = [mask mask ones(len,1)];
bwr = [blu;flipud(red)];
