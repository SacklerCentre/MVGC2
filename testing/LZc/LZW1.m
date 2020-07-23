function [n,dict] = LZW1(s)

assert(ischar(s),'input must be a character array');
s = s(:)'; % enforce row character array

dict = {}; % empty cell array
w = '';    % empty word
for c = s  % iterate through characters
	wc = [w c];
	if any(strcmp(dict,wc)) % is word in dictionary?
		w = wc;
	else
		dict = [dict;wc];   % add word to dictionary
		w = c;
	end
end

n = length(dict);
