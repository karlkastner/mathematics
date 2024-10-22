% 2024-01-24 13:56:40.759765634 +0100
% Karl Kastner, Berlin
function test

javaaddpath('./');
addpath('./mg');
dt = 1;
o = 2/3;
rng(0);

nvar = 1;
n = 16*[1,1];
L = n;
dx = L./n;

a = 2;
ad = zeros(2,3);
d = 1*[1;2];

xi = zeros(n);
x  = rand(n);
b  = rand(n);

mg_m = Multigrid();
mg_m.init({a},dt*ad,dt*d,L,n,nvar);
mg_m.o = o;
mg_m.s(1).x = x;
mg_m.s(1).b = b;
mg_m.cycle(1); %b,x,dx);
x_m = mg_m.s(1).x;

mg_2 = Multigrid();
mg_2.init_fun(@fun,L,n,nvar,xi);
mg_2.o = o;
mg_2.s(1).x = x;
mg_2.s(1).b = b;
mg_2.cycle(1); %b,x,dx);
x_2 = mg_m.s(1).x;

if (0)
mg_j = javaObject('MG');
mg_j.init(a,L,n);
mg_j.o = o;

%x = 0.*rand(n);
%b = 0.*rand(n);
%b(end/2+1,end/2+1)=1;


x_j = mg_j.cycle_(b(:), x(:), n, 0);
x_j = reshape(x_j,n);
end
x_m
x_2
x_2 - x_m

%x_j
%x_j - x_m

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

