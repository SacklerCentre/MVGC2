function [S,f] = me_pwelch(X,fs,concat,window,noverlap,nfft)

% Calculate multi-epoch autopower spectrum using Welch method

if nargin < 6, nfft     = [];    end
if nargin < 5, noverlap = [];    end
if nargin < 4, window   = [];    end
if nargin < 3, concat   = false; end % concatenate epochs?

[~,nobs,nepochs] = size(X);
if concat
	if isempty(window), window = round(nobs/4); end;
	[S,f] = pwelch(X(:,:)',window,noverlap,nfft,fs);
else
	for i = 1:nepochs
		[SS(:,:,i),f] = pwelch(X(:,:,i)',window,noverlap,nfft,fs);
	end
	S = mean(SS,3);
end
