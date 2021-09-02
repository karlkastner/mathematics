x = innerspace(0,1,1e3)'; f=20; y = sin(2*pi*x*f); F=fft(eye(length(y))); [v,e] = eig(F); y=(v*diag(diag(e).^0.5)*inv(v))*y; plot([real(y),imag(y)])

