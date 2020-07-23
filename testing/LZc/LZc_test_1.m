% mex -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Wall -Werror -O3 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -static" C/LZW_mex.c -lglib-2.0 -outdir mex

%s = char(96+randi(d,1,m));
s = char(48+(rand(1,m)>d));
%s = char(48*ones(1,m));

tic
[n,dict] = LZc(s,true);
toc

%sort(dict)
n

wlen = zeros(n,1);
for i = 1:n
	wlen(i) = length(dict{i});
end
maxlen = max(wlen)

wmaxlen = ceil(sqrt(2*length(s)))

alph = unique(s);
k = length(alph);
p = zeros(k,1);
for i = 1:k
	p(i) = nnz(s == alph(i));
end
p = p/sum(p);
entro = -nansum(p.*log(p))
