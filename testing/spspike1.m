% Locate a spectral power spike within a range of frequencies
%
% Returns NaN on failure

function ffit = spspike1(x,fs,frange,ftol)

assert(isvector(x) && isnumeric(x),'Time series must be a numeric vector');

if isempty(fs),                  fs    = 2*pi;      end % If no sampling frequency supplied assume angular frequency.
if nargin < 5 || isempty(ftol),  ftol  = sqrt(eps); end % Default tolerance for 'fminbnd' frequency fit

T = length(x);      % time series length

x = x(:)-mean(x);   % convert to column vector and temporal de-mean
t = (0:T-1)';       % time sequence

wrange = 2*pi*frange/fs; % convert to angular frequency

% Find a rough global MSE minimum, based on

m = round(T/100);

W = linspace(wrange(1),wrange(2),m);
E = -(mean(x.*cos(t*W)).^2+mean(x.*sin(t*W)).^2);
[~,i] = min(E);
wtarg = W(i);

ffit1 = fs*wtarg/(2*pi); % convert result back to ordinary frequency

w1 = wtarg-1/m;
w2 = wtarg+1/m;

% Find optimal fit frequency in range using 'fminbnd'

[wfit,~,exitflag,info] = fminbnd(@MSE,w1,w2,optimset('TolX',ftol)); % find w which minimises the MSE
assert(exitflag > 0,'Error in ''fminbnd''');
if exitflag ~= 1
	ffit = NaN;
	return
end

info

ffit = fs*wfit/(2*pi); % convert result back to ordinary frequency

ffit1
fprintf('%g\n',abs(50-ffit1))
ffit
fprintf('%g\n',abs(50-ffit))

%% Nested function to calculate MSE
%
% Passed as parameter to 'fminbnd' (see above). The MSE is calculated up to
% additive/multiplicative factors which don't depend on w.

    function mse = MSE(W)

		mse = -(mean(x.*cos(W*t))^2+mean(x.*sin(W*t))^2);

    end % function MSE

end % function sdetrend
