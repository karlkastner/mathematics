% 2022-03-13 13:13:56.621685379 +0100

n = 10;
H = [1,0.1,0.01];
%z = 2*max(H) + rand(3*n,1);
z = 2 + rand(3*n,1);
d={};
rmsd = [];
for jdx=1:length(H)
h = H(jdx);
% avoid division by zero
t = 0;
p = struct();
p.n=10;
p.L = 1;
%p.k = 5;
r = Rietkerk(p);
r.init();


for idx=1:length(z);
	z_ = z;
	z_(idx) = z_(idx) - h;
	zl = r.dz_dt(0,z_);
	z_ = z;
	z_(idx) = z_(idx) + h;
	zr = r.dz_dt(0,z_);
	
	A(:,idx) = (zr-zl)./(2*h);
end

A_ = r.jacobian(t,z);

figure(1)
clf
subplot(2,2,1)
imagesc(A)
subplot(2,2,2)
imagesc(A_)
subplot(2,2,3)
imagesc((A-A_))
colorbar
subplot(2,2,4)
imagesc(log(abs(A-A_)))
colorbar
rmsd(jdx) = rms(flat(A-A_))
d{jdx} = A-A_;
end
figure(2)
clf
subplot(2,2,1)
loglog(H,rmsd)
subplot(2,2,2)
imagesc(abs(d{end})-abs(d{end-1}))
colorbar
