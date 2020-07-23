function s = pwnoise(m,fs,mtheta,wtheta,gain)

t = (0:m-1)';

phi = mvfilter([],1,wtheta*randn(1,m))';
s = sin(2*pi*mtheta*(t-phi)/fs);

%theta = mtheta*(1 + wtheta*randn(m,1));
%s = sin(2*pi*theta.*t/fs);
