% Locate a spectral power spike within a range of frequencies

function ffit = spspike(x,fs,frange,ftol)

assert(isvector(x) && isnumeric(x),'Time series must be a numeric vector');

if isempty(fs), fs = 2*pi; end % If no sampling frequency supplied assume angular frequency.
if nargin < 4 || isempty(ftol),  ftol  = sqrt(eps); end % Default tolerance for 'fminbnd' frequency fit

T = length(x);      % time series length

x = x(:)-mean(x);   % convert to column vector and temporal de-mean
t = (0:T-1)';       % time sequence

wrange = 2*pi*frange/fs; % convert to angular frequency

% Find optimal fit frequency in range using 'fminbnd'

[wfit,~,exitflag,info] = fminbnd(@MSE,wrange(1),wrange(2),optimset('TolX',ftol)); % find w which minimises the MSE
assert(exitflag > 0,'Error in ''fminbnd''');
if exitflag ~= 1
	ffit = NaN;
	return
end

ffit = fs*wfit/(2*pi); % convert result back to ordinary frequency

%% Nested function to calculate MSE
%
% Passed as parameter to 'fminbnd' (see above). The MSE is calculated up to
% additive/multiplicative factors which don't depend on w.

    function mse = MSE(W)

		mse = -(mean(x.*cos(W*t))^2+mean(x.*sin(W*t))^2);

    end % function MSE

end
