clf
 m=1e2;
 n = 1e2;
 x = randn(n,1);
 s = abs(fft(x)).^2;
 s=s(1:n/2);
 s=s/n; %/sum(s);
 p = (1:n/2)'/(n/2+1);
 subplot(2,2,1)
 plot(1/2*chi2inv(p,2),sort(s));
 x = randn(n,m);
 s = abs(fft(x)).^2;
 s=s(1:n/2,:);
% s = s./n; %sum(s);
 b = s./mean(s,2);
 b=b(:);
 n_=numel(b);
 p=(1:n_)'/(n_+1);
 subplot(2,2,2)
 plot(sort(b),[1/2*chi2inv(p,2),m*betainv(p,1,m-1)])
 hold on
 plot([0,max(b)],[0,max(b)])
if (1)
x = randn(n,m);

 s = abs(fft(x,n*m)).^2;
 s=s(1:n*m/2,:);
 s = s./sum(s);
 b = s./mean(s,2);
 b=b(:);
 n=length(b);
 p=(1:n)'/(n+1);
subplot(2,2,3)
 plot(sort(b),[1/2*chi2inv(p,2),m*betainv(p,1,m-1)])
 hold on
 plot([0,max(b)],[0,max(b)])
end
