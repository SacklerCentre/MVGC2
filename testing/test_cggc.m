group = {[4 3] [9 8] [6 7 5 1 2]};


%-------------------------------------------------------------

g = length(group);

VARA = var9_test;
[n,~,p] = size(VARA);
V = corr_rand(n,1);

[A,C,K,info] = var_to_ss(VARA,V);
fres = info.fres;
h = fres+1;

F1 = nan(g,1);
for i = 1:g
	F1(i) = ss_to_cggc(A,C,K,V,group{i});
end

F2 = nan(g,1);
for i = 1:g
	F2(i) = var_to_cggc(VARA,V,group{i});
end
F2

F3 = nan(g,1);
for i = 1:g
	u = group{i};
	nu = length(u);
	Fi = 0;
	for k = 1:nu
		x = u(k);
		y = u; y(k) = [];
		Fi = Fi+ss_to_mvgc(A,C,K,V,x,y);
	end
	F3(i) = Fi;
end
F3

fprintf('\n');
fprintf('|F1-F3| = %g\n',maxabs(F1-F3));
fprintf('|F2-F3| = %g\n',maxabs(F2-F3));

f1 = nan(g,h);
for i = 1:g
	f1(i,:) = ss_to_scggc(A,C,K,V,group{i},fres);
end
G1 = bandlimit(f1,2);

f2 = nan(g,h);
for i = 1:g
	f2(i,:) = var_to_scggc(VARA,V,group{i},fres);
end
G2 = bandlimit(f2,2);

f3 = nan(g,h);
for i = 1:g
	u = group{i};
	nu = length(u);
	fi = zeros(1,h);
	for k = 1:nu
		x = u(k);
		y = u; y(k) = [];
		fi = fi+ss_to_smvgc(A,C,K,V,x,y,fres);
	end
	f3(i,:) = fi;
end
G3 = bandlimit(f3,2);

fprintf('\n');
fprintf('|f1-f3| = %g\n',maxabs(f1-f3));
fprintf('|f2-f3| = %g\n',maxabs(f2-f3));

fprintf('\n');
fprintf('|F1-G1| = %g\n',maxabs(F1-G1));
fprintf('|F2-G2| = %g\n',maxabs(F2-G2));
fprintf('|F3-G3| = %g\n',maxabs(F3-G3));

