% Thu  4 May 12:35:41 CEST 2017

% Test even
F=ifft(eye(n));
F(:,(n)/2+(0:2))
% Test odd
n=11; F=ifft(eye(n));
F(:,(n-1)/2+(0:3)) = []
Fi=idftmtx_man(n,3)
F-Fi
