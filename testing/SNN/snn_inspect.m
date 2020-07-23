% clear; load ~/data/SNN/modcon3.mat

[X,FS] = downsample(P,ds,fs);

[n,m] = size(X);
t = (0:m-1)/FS;

m1 = floor(t1*FS); if m1 < 1, m1 = 1; t1 = t(1);   end;
m2 = ceil(t2*FS);  if m2 > m, m2 = m; t2 = t(end); end;

t = t(m1:m2);
X = X(:,m1:m2);

gpcmds = sprintf('set xlabel"seconds"\nset ylabel "LFP"\nset xr[%f:%f]',t1,t2);
gp_qplot(t',X',[],gpcmds);
