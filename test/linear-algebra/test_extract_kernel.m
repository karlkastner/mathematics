% Mon 16 Sep 16:31:05 CEST 2024
% Karl Kastner, Berlin
n = [8,4];
[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,(n-1),2,{'circular','circular'});

id = [3,3];

K=extract_kernel(Dx,[n,n],id)
K=extract_kernel(Dy,[n,n],id)
K=extract_kernel(D2x,[n,n],id)
K=extract_kernel(Dxy,[n,n],id)
K=extract_kernel(D2y,[n,n],id)
