path(path,'../')
 n=1e4; x=pi*sort(rand(n,1)); y=sin(x(2:end-1)); L =laplacian_non_uniform(x); z=L*y;
 w = x(2:end-1).*z;
 m = z;
clf
subplot(2,1,1)
%semilogy([abs(fft(z)) abs(fft(m)) abs(fft(m2)) abs(fft(y))])
for idx=1:10
 m_ = m./x(2:end-1);
 semilogy(abs(fft(m_)),'color',rand(3,1)); hold on
 m=([w(2:end); 0] + [0; w(1:end-1)] + w)/3;
 w=m;
 %m2=([x(3:end-1).*m2(2:end); 0] + [0; x(2:end-2).*m2(1:end-1)])./x(2:end-1)/3 +m2/3;
end
legend(num2str((1:10)'))
 a=L*m;
subplot(2,1,2)
 plot(x(4:end-3),[y(3:end-2) -z(3:end-2)  a(3:end-2)]), legend('y','z','m') 

