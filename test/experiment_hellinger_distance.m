fe = fft(e); fe =fe/rms(fe); be = sqrt(S).*fe; rms(diff(be,[],2)), rms(diff(sqrt(S),[],2)), b = ifft(be); b = b.*rms(be)./rms(b); rms(diff(b,[],2))
