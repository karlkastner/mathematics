% Mon 16 Sep 16:25:47 CEST 2024
% Karl Kastner, Berlin

n = [16,16];

x = rand(n);

[R1,R1t,R12] = downsampling_matrix_2d(n);

x_ = R1*x*R1t;

I = eye(n(1)*n(2));

I=(R12*I*(4*R12'));
full(I(1:10,1:10))
K = extract_kernel(I,n/2,[3,3])

[Dx,Dy,D2x,Dxy,D2y]=derivative_matrix_2d(n,n-1,2,{'circular','circular'});
D2 = [D2x + D2y];

D2=4*(R12*D2*R12');
full(D2(1:10,1:10))
K = extract_kernel(D2,n/2,[3,3])
%x__ = reshape(R*x(:),n/2);
%x__ = reshape(R*x(:),n/2);



