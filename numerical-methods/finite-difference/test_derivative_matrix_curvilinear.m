% Mon 12 Mar 17:17:57 CET 2018

n = [70,210];
%n = [4,5];

syms x y
var = {'x','y'}
fun.z  = x.^3 + x*y.^2 + y.^4; % + x.^2 + y.^3;
fun.dz_dx   = matlabFunction(diff(fun.z,x),'var',var);
fun.dz_dy   = matlabFunction(diff(fun.z,y),'var',var);
fun.d2z_dx2 = matlabFunction(diff(fun.z,x,2),'var',var);
fun.d2z_dxy = matlabFunction(diff(diff(fun.z,x),y),'var',var);
fun.d2z_dy2 = matlabFunction(diff(fun.z,y,2),'var',var);
fun.z = matlabFunction(fun.z,'var',var);
if(1)
L = [3,3];
r     = linspace(0.5,1.5,n(1));
theta = linspace(0,2*pi,n(2));

x = cvec(r)*rvec(cos(theta));
y = cvec(r)*rvec(sin(theta));
else
	L = [3,3];
	[x,y] = meshgrid(linspace(-L(1)/2,L(1)/2,n(2)),linspace(-1.5,1.5,n(1)));
end
z        = fun.z(x,y);
dz_dx{1} = fun.dz_dx(x,y);
dz_dy{1} = fun.dz_dy(x,y);
d2z_dx2{1} = fun.d2z_dx2(x,y);
d2z_dy2{1} = fun.d2z_dy2(x,y);
d2z_dxy{1} = fun.d2z_dxy(x,y);

%Lz{1}     = d2z_dx2{1} + d2z_dy2{1};
%z        = fun.z(x,y);
%dz_dx{1} = 2.*x.^1; % + 2*x;
%dz_dy{1} = 3.*y.^2; %3*y.^2;
%d2z_dx2{1} = 2*y.^0;
%d2z_dy2{1} = 6*y;
%max(max(dz_dy{1}))
%Lz{1}    = 2 + 6*y;
m(2)     = max(flat(dz_dy{1}));

isorthogonal = false;
%[Dx,Dy,D2x,Dxy,D2y,L] = derivative_matrix_curvilinear(x,y,isorthogonal);
[Ds,Dn,D2s,Dsn,D2n] = derivative_matrix_2d(n,L);
[Dx,Dy,D2x,Dxy,D2y,L]  = derivative_matrix_curvilinear_2(x,y,isorthogonal);
spnorm(Ds(:)-Dy(:))
spnorm(Dn(:)-Dx(:))
spnorm(D2s(:)-D2y(:))
spnorm(D2n(:)-D2x(:))
spnorm(flat(Dsn-Dxy))

%Dxy = Dsn;
%full(Ds)
%full(Dx)
%pause
%Dxy=Dx*Dy;
%D2x=Dx*Dx;
%D2y=Dy*Dy;
%L = D2x+D2y;
dz_dx{2} = reshape(Dx*flat(z),n);
dz_dy{2} = reshape(Dy*flat(z),n);
%Lz{2}    = reshape(L*flat(z),n);
d2z_dx2{2} = reshape(D2x*flat(z),n)
d2z_dxy{2} = reshape(Dxy*flat(z),n)
d2z_dy2{2} = reshape(D2y*flat(z),n)
Lz{2}    = real(reshape(L*flat(z),n));

% TODO convergence test
figure(1);
clf()
figure(2);
clf()
for idx=1:2
	figure(1);
	subplot(3,3,1+3*(idx-1));
	%imagesc(dz_dx{idx})
	%plot3(x,y,dz_dx{idx})
	scatter3(flat(x),flat(y),flat(dz_dx{idx}),[],flat(dz_dx{idx}),'.')
	view(0,90)
	caxis([-1,1]*max(flat(dz_dx{1})));
	title('dz_dx')
	
	subplot(3,3,2+3*(idx-1));
	%imagesc(dz_dy{idx})
	scatter3(flat(x),flat(y),flat(dz_dy{idx}),[],flat(dz_dy{idx}),'.')
	view(0,90)
	title('dz_dy')
	caxis([-1,1]*max(flat(dz_dy{1})));

	figure(2);
	subplot(3,3,1+3*(idx-1));
	%imagesc(dz_dx{idx})
	%plot3(x,y,dz_dx{idx})
	scatter3(flat(x),flat(y),flat(d2z_dx2{idx}),[],flat(d2z_dx2{idx}),'.')
	view(0,90)
	caxis([-1,1]*max(flat(d2z_dx2{1})));
	title('d2z_dx2')
	
	subplot(3,3,2+3*(idx-1));
	%imagesc(dz_dy{idx})
	scatter3(flat(x),flat(y),flat(d2z_dy2{idx}),[],flat(d2z_dy2{idx}),'.')
	view(0,90)
	title('d2z_dy2')
	caxis([-1,1]*max(flat(d2z_dy2{1})));

	subplot(3,3,3+3*(idx-1));
	%imagesc(Lz{idx})
	scatter3(flat(x),flat(y),flat(d2z_dxy{idx}),[],flat(d2z_dxy{idx}),'.')
	caxis(limits(flat(d2z_dxy{1})));
	view(0,90)
	title('d2z_dxy')
if(0)
	subplot(3,3,3+3*(idx-1));
	%imagesc(Lz{idx})
	scatter3(flat(x),flat(y),flat(Lz{idx}),[],flat(Lz{idx}),'.')
	caxis(limits(flat(Lz{1})));
	view(0,90)
	title('laplacian')
end

end
figure(1)
subplot(3,3,7)
d= dz_dx{1}-dz_dx{2};
imagesc(d);
rms(d(:))

subplot(3,3,8)
d=dz_dy{1}-dz_dy{2};
imagesc(d);
rms(d(:))

figure(2)
subplot(3,3,7)
d= d2z_dx2{1}-d2z_dx2{2};
imagesc(d);
rms(d(:))
subplot(3,3,8)
d=d2z_dy2{1}-d2z_dy2{2};
imagesc(d);
rms(d(:))

subplot(3,3,9)
d= d2z_dxy{1}-d2z_dxy{2};
imagesc(d);
rms(d(:))
d = d(2:end-1,2:end-1);
rms(d(:))
if (0)
subplot(3,3,9)
d=Lz{1}-Lz{2};
imagesc(d);
rms(d(:))
d = d(2:end-1,2:end-1);
rms(d(:))
end


