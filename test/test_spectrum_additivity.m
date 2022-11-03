x = rand(10,1); y = rand(10,1e6); y=y-mean(y); [mean(abs(fft(x+y)).^2,2) -  (mean(abs(fft(x)).^2 + abs(fft(y)).^2 + 0*2*real(fft(y).*conj(fft(x))) + 0*conj(fft(y)).*fft(x),2))]
x = (1:10)'; y = rand(10,1e6); y=y-mean(y); 
[mean(abs(fft(x+y)).^2,2) -  (mean(abs(fft(x)).^2 + abs(fft(y)).^2 + 0*2*real(fft(y).*conj(fft(x))) + 0*conj(fft(y)).*fft(x),2))]
