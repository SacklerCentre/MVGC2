function y = gendiff_ver(x,d,r)

% Generalised differencing of vector. Calculates:
%
%     y(t) = (1-L)^(-d) x(t)
%
% where B is the lag (backshift) operator. The differencing
% order d may be fractional.
%
% INPUTS
%
% x      input vector
% d      differencing order (can be fractional)
% r      number of MA coefficients to calculate
%
% OUTPUTS
%
% y      differenced vector

if nargin < 4 || isempty(lgam), lgam = false; end

% Generalised (fractional) differencing

assert(iscolumn(x),'x must be a column vector')

% Calculate the MA coefficients


b = zeros(r+1,1);
dd = d-1;
b(1) = 1;
for k = 1:r
	b(k+1) = (1+dd/k)*b(k);
end

% Apply MA filter

y = conv(b,x,'full'); % full convolution
y = y(1:length(x));   % truncate incomplete tail (but keep incomplete head!)
