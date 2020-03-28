% Thu Aug  9 16:36:04 MSK 2012
% Karl Kästner, Berlin

syms x y z
X = [ x y z ];
%X = [ 2 3 5;
%      7 11 13;
%     17 19 23; ];
%X = rand(10,3)

p=4;
n = (p+3)*(p+2)*(p+1)/6
m = (p+2)*(p+1)*(p)/6

V = vander_3d(X,p).'
% invert ...
pause
[C_dx C_dy C_dz] = derivative_3d(V,p)


V_ = zeros(size(X,1),n);
V_ = FEM.vander_3d(V_,X,p)
V__ = zeros(size(X,1),n);
V__ = FEM.vander_3d_(V__,X,p)
%sort(V_,2)  - sort(V__,2)




C_dx_ = zeros(size(X,1),m);
C_dy_ = zeros(size(X,1),m);
C_dz_ = zeros(size(X,1),m);
C_dx_ = FEM.derivative_3d(V_, C_dx_, C_dy_, C_dz_, p)
C_dx__ = zeros(m,size(X,1));
C_dy__ = zeros(m,size(X,1));
C_dz__ = zeros(m,size(X,1));
%V__ = V__
C_dx__ = FEM.derivative_3d_(V__', C_dx__, C_dy__, C_dz__, 2)

sort(C_dx_',2)  - sort(C_dx__,2)

% nv = nv-1

