function [S,f] = me_pmtm(X,fs,concat,nw,nfft)

% Calculate multi-epoch autopower spectrum using multitaper method

if nargin < 5, nfft   = [];    end
if nargin < 4, nw     = [];    end
if nargin < 3, concat = false; end % concatenate epochs?

[~,nobs,nepochs] = size(X);
if concat
	[S,f] = pmtm(X(:,:)',nw,nfft,fs);
else
	for i = 1:nepochs
		[SS(:,:,i),f] = pmtm(X(:,:,i)',nw,nfft,fs);
	end
	S = mean(SS,3);
end
