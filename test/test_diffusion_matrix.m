 r = 10;
 a=deg2rad(-45/2);
 n=200*[1,1];
 L=n-1;
 D2 = diffusion_matrix_2d_anisotropic2(n,L,a);
  x = randn(size(D2,1),1);
 I = speye(size(D2));
 y=(I - r*D2) \ x;
 y = reshape(y,n);
 f=fftshift(abs(fft2(y)).^2);
 imagesc(f);
 axis equal

