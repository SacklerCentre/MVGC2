group = {[4 3] [9 8] [6 7 5 1 2]};

%-------------------------------------------------------------

g = length(group);


VARA = var9_test;
[n,~,p] = size(VARA);
V = corr_rand(n,1);

[A,C,K,info] = var_to_ss(VARA,V);
fres = info.fres;

F1 = ss_to_gwcgc(A,C,K,V,group);
F1

F2 = var_to_gwcgc(VARA,V,group);
F2

F3 = nan(g);
for i = 1:g
	for j = 1:g
		if j == i, continue; end
		F3(i,j) = ss_to_mvgc(A,C,K,V,group{i},group{j});
	end
end
F3

fprintf('\n');
fprintf('|F1-F3| = %g\n',maxabs(F1-F3));
fprintf('|F2-F3| = %g\n',maxabs(F2-F3));

f1 = ss_to_sgwcgc(A,C,K,V,group,fres);
G1 = bandlimit(f1,3);

f2 = var_to_sgwcgc(VARA,V,group,fres);
G2 = bandlimit(f2,3);

f3 = nan(g,g,fres+1);
for i = 1:g
	for j = 1:g
		if j == i, continue; end
		f3(i,j,:) = ss_to_smvgc(A,C,K,V,group{i},group{j},fres);
	end
end
G3 = bandlimit(f3,3);

fprintf('\n');
fprintf('|f1-f3| = %g\n',maxabs(f1-f3));
fprintf('|f2-f3| = %g\n',maxabs(f2-f3));

fprintf('\n');
fprintf('|F1-G1| = %g\n',maxabs(F1-G1));
fprintf('|F2-G2| = %g\n',maxabs(F2-G2));
fprintf('|F3-G3| = %g\n',maxabs(F3-G3));
