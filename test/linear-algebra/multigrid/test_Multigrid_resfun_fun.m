% 2024-01-24 11:48:42.916893349 +0100
% Karl Kastner, Berlin
function test

javaaddpath('./');
addpath('./mg');
dt = 100;
rng(0)

nvar = 1;
n = 8*[1,1];
L = n;
dx = L./n;

a = 2;
ad = zeros(2,3);
d = 1*[1;2];

%ad = rand(2,3);
xi = zeros(n);
x = rand(n);
b = rand(n);

%x = ones(n);
%b = ones(n);

mg_m = Multigrid();
mg_m.init({a},dt*ad,dt*d,L,n,nvar);
mg_m.s(1).x = x;
mg_m.s(1).b = b;
mg_m.resfun(1);
res_m = mg_m.s(1).res;

mg_2 = Multigrid();
mg_2.init_fun(@fun,L,n,nvar,xi);
mg_2.s(1).x = x;
mg_2.s(1).b = b;
mg_2.resfun(1);
res_2 = mg_2.s(1).res;


if (0)
mg_j = javaObject('Multigrid_java');
mg_j.init(a,ad,d,L,n);
%mg_j.s(1).x = x(:);
%mg_j.s(1).b = b(:);
%res_j = zeros(n);
%mg_j.resfun(res_j(:),b(:),x(:),n);
res_j =  mg_j.resfun_(b(:)', x(:)', 0);
%res_j = mg_j.resfun(0); %b(:),x(:),n);

%res_j = 
res_j = reshape(res_j,n);
end

res_m
res_2
res_2 - res_m
rms(res_2-res_m,'all')
%res_j
%res_j - res_m

%rms(res_j-res_m)
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
	di = ones(n(1)*n(2),1)*[-dt*d(1)/dx(1)^2, ...
			        -dt*d(1)/dx(1)^2, ...
			        -dt*d(2)/dx(2)^2, ...
			        -dt*d(2)/dx(2)^2, ...
			        +dt*2*(d(1)/dx(1)^2+d(2)/dx(2)^2)+a];
	di = reshape(di,n(1),n(2),5,nvar);	
end

end
