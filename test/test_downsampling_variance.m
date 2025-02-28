% Tue  7 Nov 20:24:07 CET 2023

-> the down and upsampling does not preserve variance
-> can be fixed to [1/6,2/3,1/6] -> but upsampling is still [1/2,1,1/2]
-> is 2x2 resampling really flawed?
-> 3 : 1 is better ? (averaging of 3 points), similar to 2x2, but it preserves location
-> 

m=20;
 x = randn(2.^m,1);
 s=[];
 n=[];
 for idx=1:m;
 n(idx,1)=length(x);
 s(idx,1)=std(x);
 x = downsampling_matrix(length(x))*x;
 end;
 loglog(n,[s,sqrt(n)/sqrt(n(1))*s(1)],'.');
 w=[1/6,2/3,1/6];
 w=w./sum(w);
 s1 = sqrt(sum(w.^2));
 s(:,2) = s(1)*s1.^(0:m-1);
 loglog(n,[s,1.01*sqrt(n)/sqrt(n(1))*s(1)],'.');

fourier_upsample

