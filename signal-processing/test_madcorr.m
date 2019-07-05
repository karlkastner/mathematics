% 2016-08-03 18:35:27.045138738 +0200
function rc

m = 10;
n = 1e5;
R = [];
R = (1:m)'/(m+1);

for idx=1:m
%	r = randn(n,1);
	[x y] =randar1_dual(10,0.1,0,0,R(idx,1),n,1);
%	R(idx,2) = sqrt(median((r-median(r)).*(s-median(s))/(median(abs(s-median(s)))*median(abs(r-median(r))))));
%	R(idx,1) = R(idx,1);
	r = madcorr(x,y);
	R(idx,2) = r;
%	R(idx,3) = sqrt(0.5*(1-R(idx,2).^2));
%	R(idx,3) = corr(r,s);
%	R(idx,5) = sqrt((1+R(idx,1))/2) - sqrt((1-R(idx,1))/2);
%	R(idx,2) = std(r);
%	R(idx,3) = mad2sd(mad(r,1));
%	R(idx,3) = sqrt(1-(1-2*r^2)^2);
%	R(idx,3) = sqrt(1-(1-r^2).^2);
%	R(idx,3) = r*sqrt(2-r^2);
end
figure(1);
clf
plot(R)

end


