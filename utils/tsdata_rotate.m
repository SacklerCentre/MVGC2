function [Y,o] = tsdata_rotate(X,o)

% Left-rotate time series X by o observations. All variables
% in X are rotated by the *same* offset.
%
% If o < 0, -o is actual offset.
% If o > 0, offset is uniform random on o:(n-o)
%
% For GC permutation testing, ensure offset large enough that
% autocorrelation has decayed below statistical significance.
% Autocorrelation decays roughly exponentially, with decay rate
% given by spectral radius; see function 'decorrlags'.

assert(ismatrix(X),'Not multi-trial; time-series data must be a 2D (variables x observations) matrix');

n = size(X,2);

if o < 0 % o is absolute offset

	o = -o;
	assert(2*o <= n && o >= 1,'Offset out of range');

else     % o is minimum for random offset

	assert(2*o <= n && o >= 1,'Minimum offset out of range');
	o = randi(n-2*o+1)+o-1; % uniform in range o:n-o

end

Y = [X(:,o+1:n) X(:,1:o)]; % left-rotate
