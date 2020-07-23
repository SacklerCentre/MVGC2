function n = LZW_alt(s)

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

dict = containers.Map;      % the dictionary

w = s(1);                   % empty current word word

for c = s(2:end)            % iterate through characters
	w = [w c];              % append current character to current word
	if ~isKey(dict,w)       % if not already in dictionary...
		dict(w) = true;     % add current word to dictionary...
		w = c;              % and reset current word to current character
	end
end
keys(dict)'
n = length(dict);           % LZW complexity is size of dictionary
