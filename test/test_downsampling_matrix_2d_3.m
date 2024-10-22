% Mon 16 Sep 15:26:15 CEST 2024
% Karl Kastner, Berlin

rng(0);
n = 8
x=rand(n); 
% PP=kron(P,P);  

[P1,P2t,P12] = downsampling_matrix_2d([n,n]);

reshape(PP*x(:),4,4) - P1*x*P2t

