% Wed  4 Sep 10:58:11 CEST 2024

if (0)
rng(0);
% generate a pattern
L    = 100;
nx   = 10*L;
fx   = fourier_axis(L,nx);
fy   = fx;

fc   = 1;
Sxpc = 1;
Syc  = Sxpc;

[fx0,sx] = normalmirroredpdf_mode2par(fc,0.5*Sxpc);
[~,sy] = normpdf_mode2par(0,Syc);

Sx = 2*normalmirroredpdf(fx,fx0,sx);
Sy = normpdf(fy,0,sy);

Sxy = cvec(Sx)*rvec(Sy);
Txy = sqrt(Sxy);

e = randn(nx,nx);
b = ifft2(Txy.*fft2(e));
% strip imaginary part introduced by rounding error
b = real(b);
b = single(b);

figure(1);
clf

for idx=1:4
disp(idx)
bmsk = true(nx);
m = 100/fc;
switch (idx)
case {1}
	% mask nothing
case {2}
	% mask x
	bmsk(m+1:end,:) = false;
case {3}
	% mask y
	bmsk(:,m+1:end) = false;
case {4}
	% mask x and y
	bmsk(m+1:end,:) = false;
	bmsk(:,m+1:end) = false;
end

sp = Spatial_Pattern('b',b,'msk.b',bmsk,'L',L*[1,1]);
sp.analyze_grid();

subplot(2,2,1)
sp.plot('S.rot.xp.con')
hold on
subplot(2,2,2)
sp.plot('S.rot.y.con')
hold on

end
legend('none','x','y','x and y')

else

% mass experiment

rng(0);

Sxpc_con = [];
Sxpc_hat = [];
Sxp_hat = zeros(200,4);
Sy_hat = zeros(200,4);
Sxp_con = zeros(200,4);
Sy_con = zeros(200,4);
Sxp_hat2 = zeros(200,4);
Sy_hat2 = zeros(200,4);
Sxp_con2 = zeros(200,4);
Sy_con2 = zeros(200,4);

% generate a pattern
L    = 20;
nx   = 10*L;
fx   = fourier_axis(L,nx);
fy   = fx;

fc   = 1;
Sxpc = 1;
Syc  = Sxpc;

[fx0,sx] = normalmirroredpdf_mode2par(fc,0.5*Sxpc);
[~,sy] = normpdf_mode2par(0,Syc);

Sx = 2*normalmirroredpdf(fx,fx0,sx);
Sy = normpdf(fy,0,sy);

Sxy = cvec(Sx)*rvec(Sy);
Txy = sqrt(Sxy);
Txy = single(Txy);
nj = 100;
for jdx=1:nj
jdx

e = randn(nx,nx);
e = single(e);
b = ifft2(Txy.*fft2(e));
% strip imaginary part introduced by rounding error
b = real(b);

for idx=1:4
bmsk = true(nx);
m    = 10/fc;
switch (idx)
case {1}
	% mask nothing
case {2}
	% mask x
	bmsk(m+1:end,:) = false;
case {3}
	% mask y
	bmsk(:,m+1:end) = false;
case {4}
	% mask x and y
	bmsk(m+1:end,:) = false;
	bmsk(:,m+1:end) = false;
end

sp = Spatial_Pattern('b',b,'msk.b',bmsk,'L',L*[1,1]);
sp.analyze_grid();

Sxp_hat(:,idx) = Sxp_hat(:,idx) + sp.S.rot.xp.hat/nj;
Sy_hat(:,idx) = Sy_hat(:,idx) + sp.S.rot.y.hat/nj;
Sxp_con(:,idx) = Sxp_con(:,idx) + sp.S.rot.xp.con/nj;
Sy_con(:,idx) = Sy_con(:,idx) + sp.S.rot.y.con/nj;
Sxp_hat2(:,idx) = Sxp_hat2(:,idx) + sp.S.rot.xp.hat.^2/nj;
Sy_hat2(:,idx) = Sy_hat2(:,idx) + sp.S.rot.y.hat.^2/nj;
Sxp_con2(:,idx) = Sxp_con2(:,idx) + sp.S.rot.xp.con.^2/nj;
Sy_con2(:,idx) = Sy_con2(:,idx) + sp.S.rot.y.con.^2/nj;

Sxpc_con(jdx,idx) = sp.stat.Sc.xp.con;
Sxpc_hat(jdx,idx) = sp.stat.Sc.xp.hat;
Syc_con(jdx,idx)  = sp.stat.Sc.y.con;
Syc_hat(jdx,idx)  = sp.stat.Sc.y.hat;

end % end for idx

end % for jdx

Sxp_hat_sd = sqrt(Sxp_hat2 - Sxp_hat.^2);
Sy_hat_sd = sqrt(Sy_hat2 - Sy_hat.^2);
Sxp_con_sd = sqrt(Sxp_con2 - Sxp_con.^2);
Sy_con_sd = sqrt(Sy_con2 - Sy_con.^2);

figure(1)
clf();
for idx=1:4;
	subplot(2,4,idx);
	ksdensity(Sxpc_hat(:,idx));
	hold on;
	ksdensity(Sxpc_con(:,idx));
	subplot(2,4,4+idx);
	ksdensity(Syc_hat(:,idx));
	hold on;
	ksdensity(Syc_con(:,idx));
end
[median(Sxpc_hat); median(Sxpc_con)]
[median(Syc_hat); median(Syc_con)]

figure(2);
clf
for idx=1:4
	subplot(2,4,idx)
	plot(fx,Sxp_hat(:,idx));
	hold on
	plot(fx,Sxp_con(:,idx));
	subplot(2,4,4+idx)
	plot(fx,Sy_hat(:,idx));
	hold on
	plot(fx,Sy_con(:,idx));
end

figure(3);
clf
for idx=1:4
	subplot(2,4,idx)
	plot(fx,Sxp_hat_sd(:,idx));
	hold on
	plot(fx,Sxp_con_sd(:,idx));
	subplot(2,4,4+idx)
	plot(fx,Sy_hat_sd(:,idx));
	hold on
	plot(fx,Sy_con_sd(:,idx));
end

end
