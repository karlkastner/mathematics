% 2015-08-04 16:56:43.850252497 +0200
% Karl Kastner, Berlin

n = 2e2;
m = 1e4;
%rho=0.999;
rho = 0.975;
% randar1 seems to have a startup problem, so skipping some vales at the beginning
nskip = 0;
a=randar1(1,rho,n+nskip,m);
a=a(nskip+1:end,:);
%a=randar1(1,rho,2*n,m);
%mu
clf;
N = (0:n-1)';
%[A B]   = acf_man(A,[],[],0);
[A B]        = acfar1(rho,n,(0:n-1)',true);
[A(:,2) B]   = acfar1(rho,n,(0:n-1)',false);
[A(:,3) B]   = acfar1_2(rho,n,(0:n-1)',true);
[A(:,4) B]   = acfar1_2(rho,n,(0:n-1)',false);
%[A_ B_] = acf_man(a,n);
%A_ = acf_man2(a,n);
s2 = var(a);
mean(s2)

biased = true;
A1     = autocorr_man4(a,n-1,[],[]);
biased = false;
A2     = autocorr_man4(a,n-1,[],[],biased);

A1 = bsxfun(@times,s2,A1)/mean(s2);
A2 = bsxfun(@times,s2,A2)/mean(s2);
S1 = sqrt(wvar(repmat(s2,n-1,1)',A1')/m)';
S2 = sqrt(wvar(repmat(s2,n-1,1)',A2')/m)';
%A_(end,:) = NaN;
%A2_(end,:) = NaN;

B_ = B;
B2_ = B;

%A  = mean(A(1:n-1,:),2);
A1 = mean(A1(1:n-1,:),2);
A2 = mean(A2(1:n-1,:),2);

figure(1);
clf();
subplot(2,2,1);
%A__ = mean(B,2); A__ = A__/A__(1);
plot(N,A,'linewidth',2);
hold on;
plot(N(1:n-1),A1,'.');
plot(N(1:n-1),A2,'.');
plot(N(1:n-1),[A2+S2,A2-S2],'k-');
plot(N(1:n-1),[A1+S1,A1-S1],'k-');
'legend analytic'
'legend sample biased'
'sample unbiased'

subplot(2,2,3)
%plot([A-A1,A-A2,A1-A2]);
legend('biased sample-model','unbiased sample - model',' unbiased - biased sample')

subplot(2,2,2);
C  = cov(bsxfun(@minus,B,mean(B,2))');
c=C(:,1);
B  = mean(B(1:n-1,:),2);
B_ = mean(B_(1:n-1,:),2);
hold on
%plot(N(1:n-1),[(B+c)/B(1) B_/B_(1)],'.');
plot(N(1:n-1),[B B_],'.');

legend('e','f','e','f');
%corrected formula')

