% Tue 27 Jul 20:06:31 CEST 2021
L  = 2;
r0 = 0.5;
r0 = 40;
n  = 1000;
p = 4;
form = 'f';

s = 3;

for idx=1:3

subplot(2,2,idx);
cla
x = linspace(0,L,n)-L/2;
[fx,Tx,msk] = fourier_axis(x);
msk = true(size(msk));
S_bp = spectral_density_bp(fx,r0,L,n,p,form);
plot(fx(msk),S_bp(msk))
hold on
switch(idx)
case {1}
	n_ = s*n;
	L_ = L;
case {2}
	n_ = n;
	L_ = s*L;
case {3}
	n_ = s*n;
	L_ = s*L;
end

% increase n
x = linspace(0,L_,n_);
[fx,Tx,msk] = fourier_axis(x);
msk = true(size(msk));
S_bp = spectral_density_bp(fx,r0,L_,n_,p,form);
plot(fx(msk),S_bp(msk),'--')

end

