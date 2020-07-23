%{
clear; LZc_test_7
%}

d = 2;

s = cell(2,1);
s{1} = 'abaaabbabbaaaaababaabbbaaab';
s{2} = 'babbbabbbaaaaaaabaababab';

n = cellfun(@length,s)

s

c1 = LZc(s{1},true,d)

c2 = LZc(s{2},true,d)

[c,cmin,cmax,] = LZc(s,true,d)
