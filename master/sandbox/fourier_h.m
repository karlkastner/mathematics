%n=256;
%n=128;
%n=64;
n=511;
L0 = 10;%2*pi;
X = L0*(1:n)'/n;
%0.5*D^2 y  + 1/r y = -y

V_c = 1./X;
%V = diag(exp(-X.^2));
V = diag(sin(2*pi*X));
%Vf = fft(V);
%D2=D2(2:n,2:n);
%A=D2;
%A = 0.5*D2;%+ Vf;
%I = eye(n);
%B = fft(I);
%A
B = [];
%A = ifft(D2) + V;
%L = (n^2)/(L0^2)*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
%A = L + V;
%A = D2 + fft(V);
%E = eigs(A,B,10,'SM')
%E=eig(A) %,B);
%sort(real(E))
%sort(E)

D2 = df(n,2,0)/(L0*2*pi)^2;
dv = ifft(D2.*fft(diag(V)));
%dv = diag(ifft(D2*fft(V))); %*L0/((n-0)^2*4*pi^2));
A = diag(ifft(D2 + fft(V_c))); % not really a matrix to compute eigenvalues from
V = diag(V);
plot([V dv])
%[1./((V./dv).*(2*pi).^(2:n+1)')]
e = [V + dv dv./V];
e(1:10,:)

