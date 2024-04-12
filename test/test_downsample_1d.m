% 2024-01-04 19:45:33.994680668 +0100

n = 10;

x = rand(10,1);
D = downsampling_matrix(n,'multigrid');

xr = D*x;

xr_ = downsample_1d(x);

rms(xr-xr_)

