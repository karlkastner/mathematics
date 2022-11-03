% Wed  6 Oct 15:10:11 CEST 2021
mm = cvec(10*2.^(0:4));
nn = mm.^2;
%mm = 2;
%nn = 1e6; %mm^2;
l = 10;
DD = zeros(length(mm),1);
for idx=1:length(mm)
 for jdx=1:l
 m = mm(idx);
 n = nn(idx);
% n=1e4;
% m=sqrt(n);
 y = randn(n,1);
% activate this to test for making the pattern instationary
 y(1:n/2) = meanfilt1(y(1:n/2),100);
 [p_(idx,1),D,pp,ratio,mdx] = periodogram_test_stationarity(y,m);
 DD(idx,1) = DD(idx,1) + D/l;
 end
end
'rmsd r'
rms(diff(ratio,[],2))
'rmsd pp'
rms(diff(pp,[],2))

subplot(2,2,1);
plot(pp(:,1),ratio);
subplot(2,2,2);
plot(pp);
vline(mdx);

p=0.01;
kolmcdf(1-p)
subplot(2,2,3)
plot(DD.*cvec(mm))

disp('D = mad(P)')
disp(DD)
disp('sqrt(n)*D')
disp(sqrt(nn).*DD)
disp('p_ = 1 - kolmcdf sqrt(n) D')
%p_ = 1 - kolmcdf(sqrt(nn).*DD);
disp(p_)
p=[0.01,0.05,0.1]
disp('kolminv(1-p)');
%disp(kolminv(p))
sqrtnD=kolminv(1-p);
disp(sqrtnD)
kolmcdf(sqrtnD)

rms(ratio(:,1)-ratio(:,2))

