% mex -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Wall -Werror -O3 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -static" C/LZW_mex.c -lglib-2.0 -outdir mex

%s = char(96+randi(d,1,m));
s = char(48+(rand(1,m)>d));
%s = char(48*ones(1,m));

tic
n = LZc(s,false)
elt1 = toc

tic
n = LZc(s,true)
elt2 = toc

elt1/elt2
