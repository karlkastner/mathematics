% 2024-09-17 14:03:15.544136296 +0200
% Karl Kastner, Berlin
n=[4,4];
 [Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,n-1,2,{'circular','circular'});
 I=eye(prod(n));
 extract_diagonals(I,n)
