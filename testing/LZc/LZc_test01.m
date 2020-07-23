% n = ???
% d = ???
% T = ???

usemex = true;

s1 = char(96+randi(d,1,n)); % high complexity
display(s1);
[c1,dict1] = LZc(s1,usemex)

s2 = char(96+ones(1,n));    % minimum complexity
display(s2);
[c2,dict2] = LZc(s2,usemex)


s3 = LZc_maxcomp(n,d);    % minimum complexity
display(s3);
[c3,dict3] = LZc(s3,usemex)
