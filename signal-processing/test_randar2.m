% 2015-08-17 17:44:42.414235480 +0200
% note: for real roots there is almost no difference at all
n=100;
 [acf A r] = acfar2(0.5,-0.4,n);
 A=[acf abs(A(1)*r(1).^(0:n-1)')];
 subplot(2,1,1);
 plot((1:n)',real(A));
 subplot(2,1,2);
 semilogy(abs(A(:,1)-A(:,2))./A(:,2))

