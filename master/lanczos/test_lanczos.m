n=20; k=3;
b = rand(n,1);
x0=b;
x0 = rand(n,1); %x0=x0/norm(x0);

A = rand(n);
A=A+A';
A = A+n*diag(diag(A));

e=zeros(k,1); e(k)=1;
% Arnoldi
[Q H] = arnoldi(A,k,x0);
[A*Q(:,1:k) - Q(:,1:k)*H(1:k,1:k) - H(k+1,k)*Q(:,k+1)*e']
% Lanczos
[Q T q_ t_ alpha beta] = lanczos(A,k,x0,1,1e-12);
A*Q - Q*T - t_*q_*e'

D=diag(beta);
Q_ = [Q q_];
T_tilda = T; T_tilda(k+1,k)=t_;

e1=zeros(k+1,1); e1(1)=1;
x = A \ b;
% solve the least squares problem by minimising the residual
y = T_tilda \ (D(1)*e1)
x_k = Q*y;
%norm(A*x_k - b)
%y = (D*T_tilda) \ (D(1)*e1)
%x_k = Q*y;
%norm(A*x_k - b)
%pause
r_true = A*x_k-b;
r_1 = Q_*(T_tilda*y) - b;
r_2 = Q_*(T_tilda*y - Q_'*b);
r_3 =    (D*T_tilda*y - Q_'*b);
r_5 =    (D*T_tilda*y - D(1)*e1);
r_6 =    (D*T_tilda*y - D(1)*1/norm(b)*e1);

'||A*x_k - b||'
norm(r_true)
'||Q_*(T_tilda*y) - b||'
norm(r_1)
'Q_*(T_tilda*y - Q_''*b) (only true if x0 is b)'
norm(r_2)
'   (T_tilda*y - D(1)*e1) (s.o.)'
norm(r_3)
'Q_*(T_tilda*y - D(1)*e1) (s.o.)'
norm(r_4)
'||D*T_tilda*y - D(1)*e1||'
norm(r_5)
norm(r_6)

[ T_tilda*y Q_'*b]
