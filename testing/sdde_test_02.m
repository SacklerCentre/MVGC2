% All times in milliseconds!

n   = 100;
T   = 10*1000;
%dt = 0.01;

amean = 5;
asdev = 1;

cmean = 9;
csdev = 2;

dmean = 100;
dsdev = 20;

pcon = 0.1;
psdf = 1;
phib = 0.2;

gpterm = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(n,1);
Y = zeros(n,1);
i = 1;
while i <= n
	X(i) = 2*rand-1;
	Y(i) = 2*rand-1;
	if hypot(X(i),Y(i)) < 1, i = i+1; end
end

a = amean+asdev*randn(n,1);

psig = (4/sqrt(pi))*pcon*psdf;

d = linspace(0,2,1000)';
p = dpdecfun(d,pcon,psig);
gp_qplot(d,p);

peff = pcon*erf(2/psig);
fprintf('\nsig     = %g\n',psig);
fprintf('min sig = %g\n',(4/sqrt(pi))*pcon);
fprintf('peff    = %g\n\n',peff);

k = 0;
for i = 1:n
	for j = 1:n
		if j == i, continue; end
		d = hypot(X(i)-X(j),Y(i)-Y(j));
		if rand < dpdecfun(d,pcon,psig)
			k = k+1;
			c(k,1) = i;
			c(k,2) = j;
			if rand < phib, csign = -1; else, csign = +1; end
			c(k,3) = csign*(cmean+csdev*randn);
			c(k,4) = d*(dmean+dsdev*randn);
		end
	end
end
K = k;

fprintf('\nedges = %d\n\n',K);

gpstem = fullfile(tempdir,[mfilename '_net']);
gp_write(gpstem,[X Y]);
gp = gp_open(gpstem,gpterm,[1 1]);
fprintf(gp,'datfile = "%s.dat"\n',gpstem);
fprintf(gp,'set size square\n');
fprintf(gp,'unset key\n');
for k = 1:K
	fprintf(gp,'set arrow from first %g,first %g to first %g,first %g\n',X(c(k,2)),Y(c(k,2)),X(c(k,1)),Y(c(k,1)));
end
fprintf(gp,'plot datfile u 1:2 pt 7\n');
gp_close(gp,gpstem,gpterm,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = ceil(T/dt)+1;
T = dt*(m-1)/1000; % seconds
t = (0:m)/1000;    % seconds

x = randn(n,m+1);

fprintf('\nsdde: '); tic
y = sdde(a,c,x,dt,true);
et = toc; fprintf('%8.6f seconds\n\n',et);

gpstem = fullfile(tempdir,[mfilename '_x']);
t = downsample(t,10);
y = downsample(y,10);
gp_qplot(t',y(1:10,:)',[],sprintf('set xr [0:%g]\nunset key',T),[],[1 Inf],[],gpstem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = dpdecfun(x,p,s)

	P = (4/sqrt(pi))*(p/s)*exp(-(x/s).^2);

end
