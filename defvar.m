function defvar(name,val)

% Assign the default value val to the variable name in the caller workspace
% if the variable name is not already defined in the caller workspace.

if evalin('caller', ['~exist(''' name ''',''var'')'])
	assignin('caller',name,val);
end
