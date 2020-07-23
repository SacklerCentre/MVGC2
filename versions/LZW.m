function [n,dict] = LZW(s,use_mex,wmaxlen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Lempl-Ziv-Welch complexity
%
% INPUT
%
% s      input character string
%
% OUTPUT
%
% n      LZW complexity
% dict   The dictionary (optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% E.g.:
%
% s = char(96+randi(5,1,1000));
% display(s);
% [n,dict] = LZW(s,true);
% display(dict);
% display(n);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(ischar(s),'Input must be a character array');


if use_mex

	if nargin < 3 || isempty(wmaxlen)
		wmaxlen = ceil(sqrt(2*length(s))); % more-or-less worst-case scenario
	else
		assert(isscalar(wmaxlen) && isnumeric(wmaxlen) && wmaxlen > 1,'Bad max. word length!');
	end

	if nargout > 1
		[n,dict] = LZW_mex(s,wmaxlen);
	else
		n = LZW_mex(s,wmaxlen);
	end

	assert(n ~= -1,'Max. word length too small!');

else
	s = s(:)';                      % enforce row vector
	dict = containers.Map;      % the dictionary
	w = s(1);                   % initialise current word to first character
	for c = s(2:end)            % iterate through remaining characters
		w = [w c];              % append current character to current word
		if ~isKey(dict,w)       % if not already in dictionary...
			dict(w) = true;     % add current word to dictionary...
			w = c;              % and reset current word to current character
		end
	end
	if length(w) > 1            % is there a word "left over"?
		if ~isKey(dict,w)       % if current word not already in dictionary...
			dict(w) = true;     % add to dictionary
		end
	end
	n = length(dict);           % LZW complexity is size of dictionary
	if nargout > 1
		dict = keys(dict)';
	end

end

% Michael S's Python code:
%{
def cpr(string):
'''
Lempel-Ziv-Welch compression of binary input string, e.g. string='0010101'. It outputs the size of the dictionary of binary words.
'''
d={}
w = ''
for c in string:
	wc = w+c
	if wc in d:
		w = wc
	else:
		d[wc] = wc
		w = c
return len(d)
%}
