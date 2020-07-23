function [S,f] = xspectrum_cmt(X,fs,win,tapers,fres,verb)

[n,m,N] = size(X);

if nargin < 3 || isempty(win),    win    = round(m/wfac); end
if nargin < 4 || isempty(tapers), tapers = [3 5];         end
if nargin < 5 || isempty(fres),   fres   = ceil(m/2);     end
if nargin < 6 || isempty(verb),   verb   = false;         end

nfft = 2*fres;
if win > nfft
	fprintf(2,'WARNING: resolution too low or window too large\n');
end
h = fres+1;
f = linspace(0,fs/2,h);

sz = size(tapers);
if sz(1) == 1 && sz(2) == 2
    tapers = dpss(win,tapers(1),tapers(2))*sqrt(fs);
else
	assert(sz(1) == win,'seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end

nwins  = floor(m/win);

S = zeros(h,n,n);
for r = 1:N
    if verb, fprintf('trial %d of %d ...',r,N); end
	for w = 1:nwins
		Xwr = X(:,1+(w-1)*win:w*win,r);
		XFT = mtfft(dtx(Xwr'),tapers,nfft,fs);
		XFT = XFT(1:h,:,:);
		for i = 1:n
			for j = 1:n
				S(:,i,j) = S(:,i,j) + mean(conj(XFT(:,:,i)).*XFT(:,:,j),2);
			end
		end
	end
    if verb, fprintf(' done\n'); end
end
f = linspace(0,fs/2,h)';
S = fs*permute(S/(nwins*N),[3 2 1]);

function XFT = mtfft(X,tapers,nfft,fs)

[m, nx] = size(X);      % size of data
[m1,nt] = size(tapers); % size of tapers
assert(m1 == m,'length of tapers is incompatible with length of data');

tapers = tapers(:,:,ones(1,nx)); % add channel indices to tapers
X      = X(:,:,ones(1,nt));      % add taper indices to data
X      = permute(X,[1 3 2]);     % reshape data to get dimensions to match those of tapers
XT     = X.*tapers;              % product of data with tapers
XFT    = fft(XT,nfft)/fs;        % fft of projected data

function y = dtx(x)

m = size(x,1);
a = zeros(m,2);
a(1:m,1) = (1:m)/m;
a(1:m,2) = 1;
y = x - a*(a\x); % Remove best fit
