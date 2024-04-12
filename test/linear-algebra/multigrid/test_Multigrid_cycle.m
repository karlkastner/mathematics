% 2024-01-24 13:56:40.759765634 +0100

javaaddpath('./');
addpath('./mg');
dt = 1e3;
o = 2/3;

n = 16*[1,1];
L = n;
dx = L./n;
a = [1,dt*1];

x = 0.*rand(n);
b = rand(n);

mg_m = Multigrid();
mg_m.init(a,L,n);
mg_m.o = o;

mg_j = javaObject('MG');
mg_j.init(a,L,n);
mg_j.o = o;

%x = 0.*rand(n);
%b = 0.*rand(n);
%b(end/2+1,end/2+1)=1;

x_m  = mg_m.cycle(b,x,dx);

x_j = mg_j.cycle_(b(:), x(:), n, 0);
x_j = reshape(x_j,n);
x_m
x_j
x_j - x_m

