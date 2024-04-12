% 2024-01-04 19:41:09.390730539 +0100
n = 10;

x = rand(10,1);
D = upsampling_matrix(2*n);

xr = D*x;

xr_ = upsample_1d(x);

rms(xr-xr_)

