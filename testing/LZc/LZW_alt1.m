function n = LZW_alt1(s)

% Lempl-Ziv-Welch complexity
%
% INPUT
%
% s      input character string
%
% OUTPUT
%
% n      LZW complexity

assert(ischar(s),'Input must be a character array');

s = s(:)';                  % enforce row array

dict = {};                  % empty dictionary (cell array)
w = '';                     % empty current word word

for c = s                   % iterate through characters
	wc = [w c];             % extended word (append current character to current word)
	if any(strcmp(dict,wc)) % is extended word already in dictionary?
		w = wc;             % yes: set current word to extended word
	else
		dict = [dict;wc];   % no:  add extended word to dictionary...
		w = c;              % ...and set current word to current character
	end
end
dict
n = length(dict);           % LZW complexity is size of dictionary
