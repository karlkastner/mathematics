% Sun 12 Jun 20:40:40 CEST 2022
n1=2^7+1;
x = (0:n1-1)';
L = n1^2;
m = 1e4;


if(1)
 a = brownian_noise_1d_interleave(n1,m);
 a = brownian_noise_1d_fft(L,[n1,m]);
%n=[1e2,1e6]; k=n(1); L=0.1;  x=linspace(0,L,n(1))';
 a=brownian_noise_1d_fft(L,[n,1e3]);
%plot([std(B')',sqrt(x)]) 
% n=[1e2,1e6]; k=n(1); L=n(1);  x=linspace(0,L,n(1))';
 [a,l,D2]=brownian_noise_1d_laplacian(L,[n1,n2]);
%plot(x,[std(B')',sqrt(x)])


subplot(2,2,1)
% plot(x,[mean(a')'])
subplot(2,2,2)
% plot(x,[std(a')',sqrt(x)])
plot(std(a')')

else
	s2=brownian_noise_1d_interleave(n1,1,true);
	plot(x,[x,s2])
end
