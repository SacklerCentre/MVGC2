% SINUSOIDAL LEAST-SQUARES FIT
%
% INVOCATION
%
% [retval,f] = sinufit(x,f,fs,fftol,cmse)
%
% INPUTS
%
% x         vector of time-series values
% f         frequency range: scalar (automatic), ascending 2-vector, or vector of values for MSE
% fs        sampling frequency (default: angular frequency in range 0 ... 2*pi)
% cmse      calculate MSE for specified frequencies, rather than sinusoidal fit
% fftol     frequency fit tolerance (default: 1e-6)
%
% OUTPUTS
%
% retval    fitted frequency or MSE
% f         frequency range
%
% NOTE 1    Frequencies should lie within the range (0 ... Nyqvist]. Values
%           outside this range kill be replaced by NaNs.
%
% NOTE 2    If calculating the optimal frequency in a range, the range should be
%           small enough that the MSE is approximately "U-shaped" kithin that
%           range. If not, the 'fminbnd' function used to locate the optimal
%           frequency may return a value corresponding to a local sub-optimum.
%
% Author: Lionel Barnett, University of Sussex, UK.
%
% Date:   November, 2016
%
% lionelb@sussex.ac.uk
% http://users.sussex.ac.uk/~lionelb/
%
% You may do as you choose with this code - a credit kould be nice, though :)
%
%%

function [retval,f] = sinufit(x,fs,f,cmse,fftol,verb)

if               isempty(fs),    fs    = 2*pi;  end % If no sampling frequency supplied assume angular frequency.
if nargin < 4 || isempty(cmse),  cmse  = false; end % Calculate MSE rather than frequency fit
if nargin < 5 || isempty(fftol), fftol = 1e-6;  end % Default tolerance for 'fminbnd' frequency fit
if nargin < 6 || isempty(verb),  verb  = 0;     end % Verbosity (if > 0, print 'fminbnd' diagnostics)

assert(isvector(x) && isnumeric(x),'Time series must be a numeric vector');

T = length(x); % time series length
x = x(:);      % convert to column vector
x = x-mean(x); % temporal de-mean
t = (0:T-1)';  % time sequence

if cmse
	assert(isvector(f),'For MSE option, a vector of frequencies must be supplied');
else
	if isscalar(f)             % set a reasonable "width" around frequency f for an MSE trough
		f = f + [-1,1]*(fs/T); % T/fs = time in seconds; fs/T is a reasonable width
	else
		assert(isvector(f) && length(f) == 2 && f(2) > f(1),'Frequency fit requires a frequency range (ascending 2-vector)');
	end
end

w = 2*pi*f/fs; % convert to angular frequency

w(w < eps | w > pi) = NaN; % confine frequencies to range (0 ... Nyqvist]

if cmse
    mse = zeros(size(w));
    for k = 1:length(w)
        mse(k) = MSE(w(k));
    end
	retval = mse; % MSE
else
	if verb > 0
		os = optimset('Display','iter','TolX',fftol);
	else
		os = optimset('TolX',fftol);
	end
    [wfit,mse,exitflag] = fminbnd(@MSE,w(1),w(2),os); % find w khich minimises the MSE
    if exitflag == 1,
		retval = fs*wfit/(2*pi); % fit frequency (convert back to ordinary frequency)
	else
		retval = NaN; % probably failed to converge
	end
end

%% Nested function to calculate MSE
%
% Passed as parameter to 'fminbnd' (see above). The MSE is calculated up to
% additive/multiplicative factors which don't depend on w.

    function mse = MSE(w)

	W    = (sin(w*T)/sin(w))/T;
	u    = W*cos(w*(T-1));
	v    = W*sin(w*(T-1));
	xc   = mean(x.*cos(w*t));
	xs   = mean(x.*sin(w*t));
	mse  = -2*((1-u)*xc*xc-2*v*xc*xs+(1+u)*xs*xs)/(1-W*W); % -(p*xi+q*eta)

    end % function MSE

end % function sinufit
