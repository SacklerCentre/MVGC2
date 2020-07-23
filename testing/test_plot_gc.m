rng_seed(1234);

F = {rand(5) rand(4)/2 rand(6); rand(3) rand(7) []};

F{1,1}(5,1) = 0;
F{1,2}(1,3) = 1;

ptitle = {...
	'Title 1',...
	'Title 2',...
	'Title 3';...
	'Title 4',...
	'Title 5',...
	''
};

Fmax = [0 1 1; 3 1 0];

%plot_gc(F,ptitle,[],Fmax,1);
%plot_gc(F,ptitle,[],Fmax,'x11');

plot_gc(F{1,1},ptitle{1,1},[],0,1);
plot_gc(F{1,1},ptitle{1,1},[],0,'x11');
