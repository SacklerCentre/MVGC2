% Supply N, Q, dt, a, fip
%

LZc_root = getenv('MATLAB_LZC');

NN = (1:N)';
QQ = (1:Q)';

figure(1); clf

rdt = sqrt(dt);
adt = 1+a*dt;

A = 1-a*dt;
V = dt;
x = varfima_to_tsdata(A,[],[],V,N,1);

Nf = 100000;
Sf = 1000000;

c    = zeros(N,Q);
crnd = zeros(N,Q);
drnd = zeros(N,Q);
fprintf('\n');
for q = 1:Q
	d = q+1; % alphabet size is q+1
	fprintf('q = %2d of %2d\n',q,Q);
	[s,qlev] = LZc_quantise(x',q);
	c(:,q) = LZc_x(s,d);
	load(fullfile(LZc_root,'data',sprintf('LZc_rand_A%02d_N%d_S%d.mat',d,Nf,Sf)));
	crnd(:,q) = cmean(1:N);
	drnd(:,q) = csdev(1:N);
	subplot(Q,1,q);
	plot(NN,x);
	for k=1:q
		line([1 N],qlev(k)*[1 1],'color','r');
	end
end

drnd  = drnd./crnd;
cnorm = c./crnd;
cmin  = repmat(LZc_cmin(NN),1,Q);

[S,w] = tsdata_to_cpsd(x,false,[],[],[],[],true,true);
S1 = V*abs(1./(1-A*exp(-1i*w))).^2; % theoretical psd
figure(2); clf
loglog(w,[S,S1]);
xlim([w(2) w(end)]);

figure(3); clf

dlegend = num2str(QQ,'%2d');
ymax = 10^ceil(log10(max(crnd(:))));

subplot(2,2,1);
loglog(NN,crnd);
ylim([1 ymax]);
title(gca,'mean LZ-complexity (random strings)');
xlabel('sequence length');
ylabel('LZc');
leg = legend(dlegend,'location','northwest');
leg.Title.Visible = 'on';
title(leg,'quantiles');
grid on

subplot(2,2,2);
loglog(NN,drnd);
title(gca,'rel. dev. LZ-complexity (random strings)');
xlabel('sequence length');
ylabel('LZc');
grid on

subplot(2,2,3);
loglog(NN,c);
ylim([1 ymax]);
title('LZ-complexity (un-normalised)');
xlabel('sequence length');
ylabel('LZc');
grid on

subplot(2,2,4);
semilogx(NN,cnorm);
ylim([0 1.5]);
title('LZ-complexity (normalised)');
xlabel('sequence length');
ylabel('LZc');
line([1 N],[1 1]);
grid on
