% 2024-01-24 12:38:37.280045952 +0100
% Karl Kastner, Berlin
function test
javaaddpath('./');
addpath('./mg');
dt = 10;
o = 2/3;
rng(0);

nvar = 1;
n = 8*[1,1];
L = n;
dx = L./n;

a = 2;
ad = zeros(2,3);
ad(1,:) = [0,-1,1];
d  = 1*[1;2];

xi = zeros(n);
x  = rand(n);
b = rand(n);

mg_m = Multigrid();
mg_m.init({a},dt*ad,dt*d,L,n,nvar);
mg_m.o = o;
mg_m.s(1).x = x;
mg_m.s(1).b = b;
mg_m.jacobi_step(1);
x_m = mg_m.s(1).x;

%ad = [ad(2,:);ad(1,:)];
%ad = -ad;
% ? why this flip?
ad = [ad(:,3),ad(:,2),ad(:,1)];

mg_2 = Multigrid();
mg_2.init_fun(@fun,L,n,nvar,xi);
mg_2.o = o;
mg_2.s(1).x = x;
mg_2.s(1).b = b;
mg_2.jacobi_step(1);
x_2 = mg_2.s(1).x;

if (0)
mg_j = javaObject('MG');
mg_j.init(a,d,L,n);


%x_m  = mg_m.jacobi_step(b,x,dx);

%res = zeros(prod(n),1);
x_j = mg_jacobi_step_(b(:)', x(:)', 0);
%x_j = mg.jacobi_step_(res, b(:), x(:), n);
x_j = reshape(x_j,n);
end

x_m
x_2
x_2-x_m
if (0)
x_j
x_j - x_m
end
if (0)
x_m = upsample1(x)
%x_j = javaMethod('MG','upsample1',flat(x'),n);
x_j = mg.upsample2(flat(x),n);
x_j = reshape(x_j,[2*n(1),n(2)])
x_j - x_m
end

% function for diagonals
function di= fun(x,n)
	dx = L./n;
	di = ones(n(1)*n(2),1)*[-dt*d(1)/dx(1)^2-dt*ad(1,1)/dx(1), ...
			        -dt*d(1)/dx(1)^2-dt*ad(1,3)/dx(1), ...
			        -dt*d(2)/dx(2)^2-dt*ad(2,1)/dx(1), ...
			        -dt*d(2)/dx(2)^2-dt*ad(2,3)/dx(1), ...
			        +dt*2*(d(1)/dx(1)^2+d(2)/dx(2)^2)-dt*ad(1,2)/dx(1)-dt*ad(2,2)/dx(2)+a];
	di = reshape(di,n(1),n(2),5,nvar);	
end

end

