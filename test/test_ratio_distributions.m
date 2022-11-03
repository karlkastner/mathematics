function test

nf = 10;
p = 0.95
n = 1e6;
 xx = 0.5*mean(randn(n,nf,2).^2,3) ;

if (0)
disp('me')
x =2*median(xx,2);
[mean(x),0] 
[std(x),0]
% note, this is not valid for very large quantiles
end

disp('mu')
x = mean(xx,2);
[mean(x),1/2*gam_mean(nf,1/nf)]
[std(x), 1/2*gam_std(nf,1/nf)]
[quantile(x,p),1/2*gaminv(p,nf,1/nf)]

% arbitrary ratio:
% x1 = gamma(a1,b1) (we have x1 = chi_2^2 ~ gam(1,1)
% x2 = gamma(a2,b2) (we have x2 = mean(chi^2) ~ 1/2*gam(nf,1/nf)
% b(a1,a2) = b2 x1/(b2 x1 + b1 x2)
%          = x1/(x1 + 1/b2 x2)	% in our case

pause()
disp('x/me, no good approximation')
d1 = 2; d2=2*(nf-1);
a=1; b = nf-1;
x = xx(:,1)./median(xx,2);
%/(chi2inv(0.5,2));
%[mean(x),std(x)]
%[std(x)/mean(x) beta_std(a,b)/beta_mean(a,b) f_std(d1,d2)/f_mean(d1,d2)]
%[quantile(x,p), betainv(p,a,b)/beta_mean(a,b), finv(p,d1,d2)/f_mean(d1,d2)]
[mean(x), nf*beta_mean(a,b)*chi2inv(0.5,2), f_mean(d1,d2)*chi2inv(0.5,2)]
[std(x), nf*beta_std(a,b)*chi2inv(0.5,2) f_std(d1,d2)*chi2inv(0.5,2)] %*chi2inv(0.5,2)^0]
[quantile(x,p) nf*betainv(p,a,b)*(chi2inv(0.5,2)), finv(p,d1,d2)*chi2inv(0.5,2)]

disp('x/mu excl')
xe = xx(:,1)./mean(xx(:,2:end),2);
d1 = 2;
d2 = 2*(nf-1);
[mean(xe), f_mean(d1,d2)]
[std(xe), f_std(d1,d2)]
q = [quantile(xe,p) finv(p,d1,d2)]
fcdf(q,d1,d2)

disp('x/mu incl')
a = 1;
b = (nf-1);
x = xx(:,1)./mean(xx,2);
[mean(x),nf*beta_mean(a,b)]
[std(x),nf*beta_std(a,b)]
q=[quantile(x,p),nf*betainv(p,a,b)]
betacdf(q/nf,a,b)

if (0)
w = (nf:-1:1); w=w/sum(w);

disp('x/mu excl, weighted')
w = (nf-1:-1:1); w=w/sum(w);
xe = xx(:,1)./sum(w.*xx(:,2:end),2);
d1 = 2;
d2 = 2*sum(w).^2/sum(w.^2);
%d2 = 2*sum(w(2:end)).^2/sum(w(2:end).^2)
[mean(xe), f_mean(d1,d2)]
[std(xe), f_std(d1,d2)]
[quantile(xe,p) finv(p,d1,d2)]
%[quantile(xe,p) finv(p,d1,d2)]
disp('x/mu incl, weighted, not working')
w = (nf:-1:1); w=w/sum(w);
x = xx(:,1)./sum(w.*xx,2);
a = 1;
nf_ = sum(w(1:end)).^2/sum(w(1:end).^2)
b = nf_-1;
[mean(x),nf_*beta_mean(a,b)]
[std(x),nf_*beta_std(a,b)]
[quantile(x,p),nf_*betainv(p,a,b)]
end
end

