function erates = var2serate(A,V,fs,fbands,fastfres)

% Calculate spectral entropy rates for VAR model for specified frequency bands
%
% A, V are the VAR parameter matrices
%
% fs is sampling rate
%
% fbands is either a column vector of frequency band boundaries, a
% 2-column matrix of frequency bands, or 'std1' for conventional
% delta, theta, alpha, gamma frequency bands as used in the neuro-
% science literature; 'std2' splits gamma into high and low gamma.
%
% fastfres specifies whether to use a faster (but potentially less
% accurate) method to estimate a good spectral resolution; the
% default (true) should be fine.
%
% Results are returned in the vector erates. In the case that fbands
% is specified as a vector of band boundaries, 'std1', or 'std2', as
% a reality check you can test whether sum(erates) == logdet(V) (if
% it's out by a couple of decimal points that's fine).

if nargin < 5 || isempty(fastfres)
	fastfres = true;
end

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
fprintf('\nFrequency bands:\n');
disp(fbands);

% Find good spectral resolution

fres = var2fres(A,V,fastfres);
fprintf('Spectral resolution = %d\n',fres);

% Calculate CPSD

S = var_to_cpsd(A,V,fres);
h = fres+1;

% Calculate log-determinants of CPSD

LDS = zeros(h,1);
for i = 1:h
	LDS(i) = logdet(S(:,:,i));
end

% Integrate across specified frequency bands

erates = zeros(nfbands,1);
for b = 1:nfbands
	erates(b) =  bandlimit(LDS,[],fs,fbands(b,:));
end
fprintf('\nEntropy rates:\n');
disp(erates);
