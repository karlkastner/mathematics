% 2024-01-24 12:38:37.280045952 +0100
javaaddpath('./');
addpath('./mg');
dt = 1e3;
o = 2/3;

n = 8*[1,1];
L = n;
x = rand(n);
dx = L./n;
a = 2;
d = 3;

mg_m = Multigrid();
mg_m.init({a},d,L,n,1);

mg_j = javaObject('MG');
mg_j.init(a,d,L,n);

x = rand(n);
b = rand(n);

mg_m.s(1).x = x;
mg_m.s(1).b = b;
mg_m.jacobi_step(1);
x_m = mg_m.s(1).x;
%x_m  = mg_m.jacobi_step(b,x,dx);

%res = zeros(prod(n),1);
x_j = mg_jjacobi_step_(b(:)', x(:)', 0);
%x_j = mg.jacobi_step_(res, b(:), x(:), n);
x_j = reshape(x_j,n);
x_m
x_j
x_j - x_m
if (0)
x_m = upsample1(x)
%x_j = javaMethod('MG','upsample1',flat(x'),n);
x_j = mg.upsample2(flat(x),n);
x_j = reshape(x_j,[2*n(1),n(2)])
x_j - x_m
end

