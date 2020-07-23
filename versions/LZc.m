function [C,Cmin,Cmax,dict] = LZc(s,normalise,d,use_mex)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Lempl-Ziv complexity
%
% INPUT
%
% s           input character string
% normalise   normalise by min/max for strings of given length and alphabet size,
%             so that C in [0,1] (must supply alphabet size, d; default = false)
% d           alphabet size
% use_mex     use C version (default: MUCH faster)
%
% OUTPUT
%
% C      LZ complexity
% dict   The dictionary - cell string (optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% E.g.:
%
%{
	d = 5;
	n = 100;
	s1 = char(96+randi(d,1,n)); % high complexity
	display(s1);
	[C1,dict1] = LZc(s1,true,d);
	display(sort(dict1));
	display(C1);
	s2 = char(96+ones(1,n));    % minimum complexity
	display(s2);
	[C2,dict2] = LZc(s2,true,d);
	display(sort(dict2));
	display(C2);
%}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(normalise)
	normalise = false;
else
	if normalise
		assert(nargin > 2 && ~isempty(d),'Must supply alphabet size for normalisation');
	end
end

if nargin < 4 || isempty(use_mex)
	use_mex = true;
end

cstr = iscellstr(s);

assert(ischar(s) || cstr,'Input must be a character string or cell string');

if use_mex

	if cstr
		n = cellfun(@length,s);
	else
		n = length(s);
	end

	if nargout > 3
		[C,dict] = LZc_mex(s,n);
	else
		C = LZc_mex(s,n);
	end

else

	dict = containers.Map;          % the dictionary
	if iscellstr(s)
		for i = 1:length(s)
			w = [];                 % initialise current word
			for c = s{i}            % iterate through input characters
				w = [w c];          % append current character to current word
				if ~isKey(dict,w)   % if current word not already in dictionary...
					dict(w) = true; % add to dictionary
					w = [];         % and reinitialise
				end
			end
		end
	else
		w = [];                     % initialise current word
		for c = s                   % iterate through input characters
			w = [w c];              % append current character to current word
			if ~isKey(dict,w)       % if current word not already in dictionary...
				dict(w) = true;     % add to dictionary
				w = [];             % and reinitialise
			end
		end
	end
	C = length(dict);               % LZ complexity is size of dictionary
	if nargout > 3
		dict = keys(dict)';         % cell string of dictionary words
	end

end

if normalise || nargout > 1
	assert(nargin > 1,'Must supply alphabet size for normalisation');
	if cstr
		fprintf(2,'WARNING: normalisation for cell strings is experimental!\n');
		n = sum(cellfun(@length,s)) % treat as concatenated string... (approximate)
	else
		n = length(s);
	end
	[Cmin,Cmax] = LZc_cminmax(d,n);
	if normalise
		C = (C-Cmin)/(Cmax-Cmin);
	end
end
