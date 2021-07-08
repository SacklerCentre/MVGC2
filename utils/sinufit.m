% SINUSOIDAL LEAST-SQUARES-FIT FREQUENCY
%
% INVOCATION
%
% [ffit,f] = sinufit(x,fs,f,fftol,verb)
%
% mse  = sinufit(x,fs,f,'mse')
%
% INPUTS
%
% x         matrix of time-series values (variables x observations)
% fs        sampling frequency (default: angular frequency in range 0 ... 2*pi)
% f         frequency range: scalar (automatic), ascending 2-vector, or vector of values for MSE
% fftol     frequency fit tolerance (default: 1e-6)
% verb      verbosity
%
% OUTPUTS
%
% ffit      fitted frequency
% mse       the MSE over the frequency range
% f         actual frequency range used
%
% NOTE      If calculating the optimal frequency in a range, the range should be
%           small enough that the MSE is approximately "U-shaped" within that
%           range. If not, the 'fminbnd' function used to locate the optimal
%           frequency may return a value corresponding to a local sub-optimum.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [retval,f] = sinufit(x,fs,f,fftol,verb)

if nargin < 2 || isempty(fs), fs = 2*pi; end % If no sampling frequency supplied assume angular frequency.

cmse = nargin > 3 && ischar(fftol) && strcmpi(fftol,'mse');

assert(ismatrix(x),'Time series must be a matrix');
[n,m] = size(x);

if cmse
	assert(nargin == 4,'Too many input arguments for ''mse'' option');
	assert(nargout < 2,'Too many output arguments for ''mse'' option');
	assert(isvector(f),'For ''mse'' option, a vector of frequencies must be supplied');
	f = f(:); % ensure column vector
else
	if nargin < 4 || isempty(fftol), fftol = 1e-6;  end % Default tolerance for 'fminbnd' frequency fit
	if nargin < 5 || isempty(verb),  verb  = false; end % Verbosity (if set, print 'fminbnd' diagnostics)

	if isvector(f)
		f = f(:) + [-1,1]*(fs/m); % ensure column vector; m/fs = time in seconds; fs/m is a reasonable width
	else
		assert(ismatrix(f) && size(f,2) == 2,'Fit frequencies must be a vector or two-column matrix');
		assert(all(f(:,2) > f(:,1)),'Fit frequencies must be ascending');
	end
end
%assert(all(f(:) > eps & f(:) <= fs/2),'Frequencies must lie in range (0 .. Nyqvist]');
nf = size(f,1);

x = x-mean(x,2); % temporal de-mean
t = 0:m-1;       % time sequence
w = 2*pi*f/fs;   % convert to angular frequency

if cmse
    mse = zeros(size(w));
    for k = 1:nf
        mse(k) = MSE(w(k));
    end
	retval = mse; % MSE
else
	if verb
		os = optimset('Display','iter','TolX',fftol);
	else
		os = optimset('TolX',fftol);
	end
	retval = nan(nf,1);
	mse    = nan(nf,1);
	for k = 1:nf
		[wfit,mse(k),exitflag] = fminbnd(@MSE,w(k,1),w(k,2),os); % find wfit which minimises the MSE
		if exitflag == 1,
			retval(k) = fs*wfit/(2*pi); % fit frequency (convert back to ordinary frequency)
		end
	end
end

% Nested function to calculate MSE
%
% Passed as parameter to 'fminbnd' (see above). The MSE is calculated up to
% additive/multiplicative factors which don't depend on w.

    function mse = MSE(ww)

	W    = (sin(ww*m)/sin(ww))/m;
	u    = W*cos(ww*(m-1));
	v    = W*sin(ww*(m-1));
	xc   = mean(x.*cos(ww*t),2);
	xs   = mean(x.*sin(ww*t),2);
	mse  = -2*sum((1-u)*xc.*xc-2*v*xc.*xs+(1+u)*xs.*xs)/(1-W*W); % -sum(p*xi+q*eta)

    end % function MSE

end % function sinufit
