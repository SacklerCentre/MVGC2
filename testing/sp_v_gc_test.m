
 a = linspace(-0.99,0.8,1000)';

 c = 2;
 %b = 0.8;

 aa = a.*a;
 bb = b.*b;
 cc = c.*c;
 ab = a.*b;
 ac = a.*c;

 P = cc.*(1+ab)+(1-bb).*(1-ab+2.*k.*ac);

 p = P./((1-aa)*(1-bb).*(1-ab));

 D = (1+bb+cc)/2;

 F = log(D+sqrt(D.*D-bb));

gp_qplot(a,p);
