n = [10,10];

x = rand(n);

[R1,R1t,R1_,R2_] = downsampling_matrix_2d(n);

x_ = R1*x*R1t;
%x__ = reshape(R*x(:),n/2);
%x__ = reshape(R*x(:),n/2);



