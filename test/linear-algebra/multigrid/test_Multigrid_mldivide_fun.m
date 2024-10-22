% 2024-01-24 13:56:33.251648360 +0100
% Karl Kastner, Berlin

function test
javaaddpath('./');
rng(0);
ntime = 2;

nvar = 1;
dt     = 1;
o      = 2/3;
gamma  = 1;
reltol = 1e-7;

n = 2^8*[1,1];
L = n;
x = rand(n);
dx = L./n;

a = 1;
d = [1;1];
ad = zeros(2,3);
ad(1,:) = -10*[0,+1,-1];

vx = -1;

%ad(1,1:3) = 0*upwind_kernel(vx);
%ad(2,2:3) = 0*[-1,1];

x0 = 0.*rand(n);
b = 0.*rand(n);
b(end/2+1,end/2+1)=1;

mg_m = Multigrid();
mg_m.init({a},dt*ad,dt*d,L,n,nvar);
mg_m.o = o;
mg_m.gamma = gamma;
mg_m.reltol = reltol;
% ? why this flip?
ad = [ad(:,3),ad(:,2),ad(:,1)];

xi = zeros(n(1),n(2),nvar);
mg_2 = Multigrid();
mg_2.init_fun(@fun,L,n,nvar,xi);
mg_2.o = o;
mg_2.gamma = gamma;
mg_2.reltol = reltol;

if (0)
mg_j = javaObject('Multigrid_java');
mg_j.init(a,dt*ad,dt*[d],L,n);
mg_j.o = o;
mg_j.nsubcycle = gamma;
mg_j.reltol = reltol;
end

x_m = b;
x_2 = b;
x_j = b(:)';

for idx=1:ntime

tic()
[x_m,resn_m]  = mg_m.mldivide(x_m,x_m);
toc()

tic();
x_2 = mg_2.mldivide(x_2, x_2);
%x_j = mg_j.mldivide(x_j, x_j);
toc()
end

% x_j = reshape(x_j,n);

if (0)
x_m
x_j
x_j - x_m
end
if (0)
nr = 1e2;
T = rand(size(x));
x=b;
tic();
for idx=1:nr
	x=ifft2(T.*fft2(x));
end
printf('ET %f\n',toc()/nr);
end

printf('rms(x) %e %e\n',rms(x_m(:)),rms(x_2(:)))
printf('Error %e\n',rms(x_2-x_m,'all'));
%printf('iter m %d j %d\n',length(resn_m),mg_j.iter)
%printf('resn %e %e\n',resn_m([end]),mg_j.resn);

if (0)

figure(1)
clf
subplot(2,3,1)
imagesc(x_m)
colorbar
subplot(2,3,2)
imagesc(x_j)
colorbar

subplot(2,3,3)
imagesc(x_j-x_m)
colorbar

end

% function for diagonals
function di= fun(x,n)
	dx = L./n;
	di = ones(n(1)*n(2),1)*[-dt*d(1)/dx(1)^2, ...
			        -dt*d(1)/dx(1)^2, ...
			        -dt*d(2)/dx(2)^2, ...
			        -dt*d(2)/dx(2)^2, ...
			       +dt*2*(d(1)/dx(1)^2+d(2)/dx(2)^2)+a];
	di = reshape(di,n(1),n(2),5,nvar);	
end

end

