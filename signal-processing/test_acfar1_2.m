% 2015-08-06 17:15:26.715013075 +0200
function test_acfar1()


n = 500;
delay = 75;
nf = 10;

%R = linspace(0.5,1.1,10);
%figure(1);
%clf();
%for idx=1:length(R)
%	rat(idx) = ratio(R(idx),n,nf,delay);
%	subplot(3,3,idx)
%	plot(acfar1(R(idx),n));
%end

m=100;
R = linspace(0.999,1,m);
z = NaN(m,1);
for idx=1:length(R)
	a = acfar1(R(idx),n);
	f = find(a < exp(-1),1,'first');
	if (~isempty(f))
		z(idx) = f;
	end
end

figure(2)
clf
plot(R,z)
%vline(1);

end

function ratio = ratio(rho,n,nf,delay)
        a = acfar1(rho,n);
        ratio = a(nf+1+delay)/a(nf+1);
end


