function [n,dict] = LZW_ver(s,method)

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

s = s(:)';                      % enforce row array

switch method

case 'carray'

	dict = {};                  % empty dictionary (cell array)

	w = s(1);                   % initialise current word to first character
	for c = s(2:end)            % iterate through remaining characters
		w = [w c];              % append current character to current word
		if ~ismember(w,dict)    % if not already in dictionary...
			dict = [dict;w];    % add current word to dictionary...
			w = c;              % and reset current word to current character
		end
	end

	n = length(dict);           % LZW complexity is size of dictionary

case 'map'

	dict = containers.Map;      % the dictionary

	w = s(1);                   % initialise current word to first character
	for c = s(2:end)            % iterate through remaining characters
		w = [w c];              % append current character to current word
		if ~isKey(dict,w)       % if not already in dictionary...
			dict(w) = true;     % add current word to dictionary...
			w = c;              % and reset current word to current character
		end
	end

	n = length(dict);           % LZW complexity is size of dictionary

	if nargout > 1
		dict = keys(dict)';
	end

case 'sll'

	[n,dict] = LZW_mex_sll(s);

case 'hash'

	[n,dict] = LZW_mex_hashset(s);

otherwise

	error('Unknown method');

end
