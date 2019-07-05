% 2017-07-20 12:02:10.857655392 +0200

if (0)
	r=linspace(0,1);
	L=[10 1e2 1e3];
	for jdx=1:length(L);
		for idx=1:length(r);
			w(idx)=last(acfar1(r(idx),L(jdx),1))/r(idx)
		end;
		plot(w);
		hold on;
		ylim([0 1]);
	end
else

%n = 10;
m = 1e3;
rho = 0.1;
k = 10;
R = (89 + (1:k))/100;
w=[];
R_ = [];
N = [10:10:100];
for jdx=1:length(N)
jdx
n =N(jdx);
for idx=1:length(R)
	rho=R(idx);
	x = randar1(1,rho,n,m);
	a = autocorr_man5(x,1);
	a = mean(a(2,:));
	R_(idx,jdx) = a;
	w(idx) = a/rho;
	p(idx,jdx) = log(R(idx))/log(a);
end
%clf
%plot(R,w)
end
clf
plot(R,p)
%hold on
%pause

end

