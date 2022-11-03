% 2021-06-26 23:08:47.946762034 +0200

n=100;
 D = derivative_matrix_2_1d(2*n+1,2*n,2,'circular');
D=full(D);
 rho=0.2;
 s=0;
 for idx=0:-1;
 s=s+(rho*D)^idx;
 end;
 s=inv(speye(2*n+1)-rho/(1-2*rho+rho^2)*D);
 s=s(n+1,n+1:n+4);
 s=s/s(1);
 [s;
 rho.^(0:3)]
