% Oct  8 17:07

%FDM analogue to example of trefethen
% Programme 8 with FDM
% Trefethen: spectral met in matlab

 n=100;%36;
 L=2*8;
 x=L*(0:n-1)/n - L/2;
 A=(n^2)/(L^2)*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
% A(1,end)=A(1,2);
% A(end,1)=A(end,end-1);
% D2 = df(n,2,0);
% A = ifft(diag(D2)/L^2);
 E=eig(-A+diag(x.^2));

E(1:10)
sort(E)

