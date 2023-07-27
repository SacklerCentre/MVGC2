function [erates,fres] = var2serate(A,C,K,V,fs,fbands,rois,fres)

% Calculate spectral entropy rates for SS model for specified
% frequency bands and ROIs.
%
% A, C, K, are the ISS parameter matrices, V the residuals covariance matrix
%
% fs is sampling rate
%
% fbands is either a column vector of frequency band boundaries, a
% 2-column matrix of frequency bands, or 'std1' for conventional
% delta, theta, alpha, gamma frequency bands as used in the neuro-
% science literature; 'std2' splits gamma into high and low gamma.
% If fbands is left empty (default), then "broadband" entropy rates
% (i.e.,across the entire spectrum) are returned.
%
% rois is a cell vector, where each cell is a vector of channel
% numbers for the channels in an ROI. It may also take two special
% values: 'allchans' treats the entire system as a single roi;
% 'perchan' returns entropy rates for each individual channel.
%
% fres specifies the spectral resolution (number of frequencies,
% evenly spaced, at which spectral power is calculated). If left
% empty (default), a minimum optimal value is calculated
% automatically. If it takes the special value 'accest', then
% a slower, but sometimes more optimal method is used for auto-
% matic calculation (the default faster method is generally okay,
% though).
%
% Results are returned in the vector erates. In the case that fbands
% is specified as a vector of band boundaries, 'std1', or 'std2', as
% a reality check you can test whether sum(erates) == logdet(V) (if
% it's out by a couple of decimal points that's fine). Broadband
% entropy rates should also integrate to the corresponding logdet(V).

% Some parameter defaults

if nargin < 6
	fbands = []; % broadband
end

if nargin < 7 || isempty(rois)
	rois = 'allchans'; % return whole-system entropy rates (i.e., all channels considered as a single ROI)
end

if nargin < 8
	fres = []; % automatic spectral resolution (fast method)
	fastfres = true;
elseif ischar(fres) && strcmpi(fres,'accest')
	fres = [];
	fastfres = false;
else
	assert(isscalar(fres) && isnumeric(fres),'Spectral resolution must empty, ''accest'', or a number');
end

broadband = isempty(fbands); % return entropy rate at every frequency (at resolution fres) from 0 - Nyqvist

if broadband
	fprintf('\nBroadband: 0-%g Hz\n',fs/2);
else

	% Sort out frequency bands specification

	if ischar(fbands)
		if     strcmpi(fbands,'std1') % single gamma band > 30 Hz
			fbands = [4;8;15;30];
		elseif strcmpi(fbands,'std2') % gamma band split into low (30 - 50 Hz) and high (> 50 Hz)
			fbands = [4;8;15;30;50];
		else
			error('Unknown frequency bands spec');
		end
	else
		assert(ismatrix(fbands) && (size(fbands,2) == 1 | size(fbands,2) == 2),'Frequency bands must be a string spec or a 1 or 2 column matrix');
	end
	if size(fbands,2) == 1 % vector of band boundaries
		assert(issorted(fbands),'Frequency boundaries must be sorted ascending');
		fb = fbands;
		nbb = size(fbands,1);
		fbands = zeros(nbb+1,2);
		fbands(1,:) = [0 fb(1)]; % all lower frequencies
		for b = 2:nbb
			fbands(b,:) = [fb(b-1) fb(b)];
		end
		fbands(nbb+1,:) = [fb(nbb) fs/2]; % up to Nyqvist frequency
	end
	nfbands = size(fbands,1);
	fprintf('\nFrequency bands (Hz):\n');
	disp(fbands);
end

% Sort out ROIs

nchans = size(V,1);
if ischar(rois) && strcmpi(rois,'allchans')    % all channels as a single ROI
	roi = {1:nchans};
	perchan = false;
	nrois = 1;
elseif ischar(rois) && strcmpi(rois,'perchan') % each channel as an ROI
	perchan = true;
	nrois = nchans;
else
	assert(isvector(rois) && iscell(rois),'ROI spec must be ''allchans'', ''perchan'', or a cell vector');
	perchan = false;
	roi = rois;
	nrois = length(rois);
end

% Find good spectral resolution

if isempty(fres) % automatic calculation
	fres = ss2fres(A,C,K,V,fastfres);
end
fprintf('Spectral resolution = %d\n',fres);

% Calculate CPSD

S = ss_to_cpsd(A,C,K,V,fres);
h = fres+1;

% Calculate log-determinants of CPSD

LDS = zeros(h,nrois);
if perchan
	for i = 1:nchans
		LDS(:,i) = log(squeeze(S(i,i,:))); % log-autopwer!
	end
else
	for r = 1:nrois
		Sr = S(roi{r},roi{r},:);
		for k = 1:h
			LDS(k,r) = logdet(Sr(:,:,k));
		end
	end
end

if broadband
	erates = LDS;
else

	% Integrate across specified frequency bands

	erates = zeros(nfbands,nrois);
	for r = 1:nrois
		for b = 1:nfbands
			erates(b,r) =  bandlimit(LDS(:,r),[],fs,fbands(b,:));
		end
	end
	fprintf('\nEntropy rates:\n');
	disp(erates);
end
