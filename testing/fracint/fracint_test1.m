

a = fracint_coeffs(d,r,true);
b = fracint_coeffs(d,r,false);

A = d*realpow(0:r,-d-1);
B = d*realpow(0:r,d-1);

gp_qplot((1:r)',abs([a(2:end)' A(2:end)' b(2:end)' B(2:end)']),{'AR','AR1','MA','MA1'},'set logs xy');
