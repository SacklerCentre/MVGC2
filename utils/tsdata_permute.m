function Y = tsdata_permute(X,pparm)

% Permute multi-trial time series for permutation testing.
%
% Cyclic permutation
% ------------------
%
% If pparm > 0, then we perform a cyclic permutation, left-rotating the
% data, with trials concatenated, by a random offset. The data is
% then reassembled into multi-trial form. The parameter mino = pparm
% specifies the minimum offset permitted; this minimum offset should,
% if possible, be set large enough that autocorrelation decays below
% statistical significance. Autocorrelation decays roughly exponentially,
% with decay rate given by spectral radius; use function the 'decorrlags'
% (with n = number of observations per trial) to obtain a reasonable
% estimate of the minimum offset to guarantee decorrelation.
%
% Block derangement
% -----------------
%
% If pparm < 0, blocks of observations of size bsize = -pparm from all
% trials are 'deranged' (permuted with none left in-place) and reassembled
% into multi-trial form.
%
% NOTE 1: block size must divide number of observations per trial exactly
% (because the developers are too lazy to deal with trailing observations),
% so you may need to trim some observations off your data!
%
% NOTE 2: blocks should be large enough that they retain a good amount of
% autocovariance structure in the data, but small enough that there are a
% decent number of distinct derangements available.

[nvars,nobs,ntrials] = size(X);

if     pparm > 0 % cyclic permutation

	mino = pparm;                      % minimum rotation offset
	tnobs = ntrials*nobs;              % total number of observations
	assert(2*mino <= tnobs && mino >= 1,'Minimum offset out of range');
	o = randi(tnobs-2*mino+1)+mino-1;  % uniform in range mino:(tnobs-mino)
	X = X(:,:);                        % concatenate trials
	Y = [X(:,o+1:tnobs) X(:,1:o)];     % left-rotate data
	Y = reshape(Y,nvars,nobs,ntrials); % de-concatenate rotated data into trials

elseif pparm < 0 % block permutation

	bsize = -pparm;                     % permutation block size
	nblocks = floor(nobs/bsize);        % number of blocks per trial
	assert(nblocks*bsize == nobs,'Sorry, block size must divide number of observations per trial exactly - trim your data if necessary!');
	tblocks = nblocks*ntrials;          % total number of blocks across all trials
	Y = reshape(X,nvars,bsize,tblocks); % stack up blocks across all trials
	Y = Y(:,:,derangement(tblocks));    % derange stacked blocks (permutate leaving no block in place)
	Y = reshape(Y,nvars,nobs,ntrials);  % unstack permuted blocks into trials

else
	error('Permutation parameter must be non-zero');
end

function d = derangement(n)

% https://uk.mathworks.com/matlabcentral/fileexchange/30189-randpermfull
% Copyright (c) 2016, Jos (10584) All rights reserved.

if n == 2                    % the trivial case of two elements
	d = [2 1];
else                         % a fast rejection scheme
	baad = true;
	while baad
		d = randperm(n);     % create a possible derangement
		baad = false;        % and assume that it is :-)
		for k = 1:n          % check for derangement by looping over the elements
			if d(k) == k     % violation found
				baad = true; % try again
				break        % no need to check further
			end
		end
	end
end
