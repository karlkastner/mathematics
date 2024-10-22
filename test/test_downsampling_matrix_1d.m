% Mon 16 Sep 15:26:15 CEST 2024
% Karl Kastner, Berlin

rng(0);
n = 16
x=rand(n,1); 
% PP=kron(P,P);  

[P] = downsampling_matrix_1d(n,'mg');

mean(x)
mean(P*x)

%reshape(PP*x(:),4,4) - P1*x*P2t

I = eye(n);
I = P*I*(2*P');
I(2,:)

D2 = derivative_matrix_2_1d(n,n-1,2,'circular');
full(D2(2,:))
D2 = P*D2*(2*P');
full(D2(2,:))
size(D2)

D1 = derivative_matrix_1_1d(n,n-1,2,'circular');
full(D1(2,:))
D1 = P*D1*(2*P');
full(D1(2,:))
size(D1)
