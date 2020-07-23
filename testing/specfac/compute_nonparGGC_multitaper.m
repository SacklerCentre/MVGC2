%==========================================================================
% Written by M.Dhamala (August 2006).
% Modified by M.F.Pagnotta (June 2018).
%--------------------------------------------------------------------------
% INPUT
% - X:          bivariate or multivariate signals in the form of 3D matrix [nTime, nTrial, nChannel]
% - fs:         data sampling rate in Hz
% - df:       a desired frequency resolution (e.g. 1) - df is an optional input,  default value is fs/size(X,1)
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - flgCond:    if 1 do conditional Granger causality, while if 0 do pairwise causality
% - NW:         time-bandwidth parameter to control the number of tapers (given by 2*NW-1)
%--------------------------------------------------------------------------
% OUTPUT
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - causality:  Granger-Geweke causality between all pairs of signals (channels) [nFreqs, nChannel, nChannel]
%               e.g., causality(:,2,1) means causality from 2 to 1, and 1 to 2 is causality(:,1,2)
%               self-causality are set to zero, i.e. causality(:,k,k) = 0 for all k
% - COH:        coherence between all pairs of signals [nFreqs, nChannel, nChannel]
% - iCOH:       imaginary-coherence between all pairs of signals [nFreqs, nChannel, nChannel]
% - power:      1-sided auto-spectra [nFreqs, nChannel]
% - S:          spectral density matrix [nChannel, nChannel, nFreqs]
%==========================================================================
% If you use nonparametricGGC_toolbox for a paper or talk please include
% the following references:
%
% M.F. Pagnotta, M. Dhamala, G. Plomp, Benchmarking nonparametric Granger
% causality: robustness against downsampling and influence of spectral
% decomposition parameters, NeuroImage. 183 (2018) 478–494.
% https://doi.org/10.1016/j.neuroimage.2018.07.046
%
% M. Dhamala, G. Rangarajan, M. Ding, Analyzing information flow in brain
% networks with nonparametric Granger causality, NeuroImage. 41 (2008)
% 354–362. https://doi.org/10.1016/j.neuroimage.2008.02.020.
%
%_____________________
% Useful references:
% [1] Pagnotta et al., NeuroImage, 2018.
% [2] Dhamala et al., NeuroImage, 2008.
% [3] Dhamala et al., Phys. Rev. Lett., 2008.
% [4] Thomson, Proc. IEEE, 1982.
%
%__________________________________________________________________________
% [License]
%
% This file is part of nonparametricGGC_toolbox.
% Copyright (©) 2018, Mattia F. Pagnotta.
%
% nonparametricGGC_toolbox is free software: you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% nonparametricGGC_toolbox is distributed in the hope that it will be
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with nonparametricGGC_toolbox.  If not, see
% <https://www.gnu.org/licenses/>.
%__________________________________________________________________________
%==========================================================================

function [S,freq] = compute_nonparGGC_multitaper(X,fs,freq,NW)

%--- Initialization --------------
[Nt,Ntr,Nc] = size(X)                                                    % Nt: timepoints, Ntr: trials, Nc: channels

if nargin < 4 || isempty(NW), NW = 4; end; % by default: # of tapers is 7

df = freq(2)-freq(1);

[S,f] = sig2mTspect_nv(X,fs,NW,df);  % Not vectorized, less memory-demanding, signals to multitapered auto&cross-spectra

S = S(:, :, 1:length(freq));      % reduce S to frequency range of interest

%S(:,:,1) = [];

%S = (fs/2)*S;                                                        % omitting zero-frequency spectra (auto and cross)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function [S,f] = sig2mTspect_nv(X,fs,nw,df)
%---------------------
% This function computes auto- & cross- spectra by using multitapers.
%
% INPUT
% - X:          multichannel data [Nt, Ntr, Nc]
% - fs:         data sampling rate in Hz
% - nw:         parameter to control the number of tapers (= 2*nw-1), good nw are 1.5, 2, 3, 4, 5, 6, or 7
% - df:       a desired frequency resolution (e.g. 1) - df is an optional input,  default value is fs/size(X,1)
% OUTPUT
% - S:          matrix of auto- and cross-spectra [Nc, Nc, nFreqs]
%---------------------

[Nt,Ntr,Nc] = size(X);        % Nt: timepoints, Ntr: trials, Nc: channels

%{
if df>fs/Nt,
       npad = 0;
       df = fs/Nt;
disp('REJIG')
else
disp('NO')
end
%}
if df<=fs/Nt,
    npad = round((fs/df-Nt)/2)                % These many zeros will be padded on each side of the data
disp('PAD')
else
	npad = 0;
    npad = round((fs/df-Nt)/2)                 % These many zeros will be padded on each side of the data
disp('NO PAD')
end
npad = 0;
f = fs*(0:fix((Nt+2*npad)/2))/(Nt+2*npad);      % upto Nyquist-f

[tapers, v] = dpss(Nt+2*npad, nw, 2*nw-1);      % Return the 2*nw-1 most band-limited discrete prolate spheroidal sequences

S = zeros(Nc,Nc,Nt+2*npad);

for itrial = 1:Ntr,
    for ii = 1:Nc,
        Xft(:,:,ii) = mtfft(squeeze(X(:,itrial,ii)),tapers,fs,npad);
    end
    for ii = 1:Nc,
        for jj = 1:Nc,
            s(ii,jj,:) = squeeze(mean(Xft(:,:,ii).*conj(Xft(:,:,jj)),2)); %averaging over tapers
        end
    end
    S = S + s;
end
S = S/Ntr;                              % Averaging over trials
S = 2*S(:,:,1:fix(end/2)+1)/fs;         % factor 2 for make one-sided spectra
S(:,:,1) = S(:,:,1)/2;                  % dc-power doesn't double for one-sided case

%==========================================================================
function xf  = mtfft(data,tapers,fs,npad)
%---------------------
% This function multiply the time series from each trial by the preselected
% number of orthogonal tapers. The products are then Fourier-transformed.
% (Used inside the function sig2mTspect_nv.m)
%
% INPUT
% - data:       time series from each trial
% - tapers:     preselected orthogonal tapers
% - fs:         sampling frequency
% - npad:       number of zeros for padding
% OUTPUT
% - xf:         Fourier-transformed of the products between each trial time series and the tapers
%---------------------

x0   = zeros(npad,size(data,2));
data = cat(1,x0,data);
data = cat(1,data,x0);
data = data(:,ones(1,size(tapers,2)));
data = data.*tapers;
xf   = fft(data,[],1);
