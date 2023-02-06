n = 11;
f = [1:n] - (n+1)/2
id = fourier_freq2ind(f,n)
f(id) = f
f_ = fourier_axis(1,n)'
rms(f-f_)

n = 10;
f = [0:n-1] - (n)/2
id = fourier_freq2ind(f,n)
f(id) = f
f_ = fourier_axis(1,n)'
rms(f-f_)
