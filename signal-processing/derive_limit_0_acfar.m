% 2017-08-13 11:00:13.067548426 +0200
syms dn r w L p;
 a = acfar1(exp(-dn/L),w/dn,p*w/dn);
 a=limit(limit(a,dn,0),L,inf);
 s=solve(a-'exp(-1)',p);

syms dn r w L p;
 a = acfar1(exp(-dn/L),w/dn,p*w/dn);
 a=limit(limit(a,dn,0),L,inf);
 solve(a)

