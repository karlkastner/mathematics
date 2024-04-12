% 2024-01-24 11:48:42.916893349 +0100
javaaddpath('./');
addpath('./mg');
dt = 1e3;
rng(0)

a = 2
d = 3*[1;2];
n = 8*[1,1];
L = n;
x = rand(n);
dx = L./n;

ad = zeros(2,3);
ad = rand(2,3);

x = rand(n);
b = rand(n);
%b = ones(n);

mg_m = Multigrid();
mg_m.init({a},ad,d,L,n,1);
mg_m.s(1).x = x;
mg_m.s(1).b = b;
mg_m.resfun(1);
res_m = mg_m.s(1).res;

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
res_m
res_j
res_j - res_m

rms(res_j-res_m)
if (0)
x_m = upsample1(x)
%x_j = javaMethod('MG','upsample1',flat(x'),n);
x_j = mg.upsample2(flat(x),n);
x_j = reshape(x_j,[2*n(1),n(2)])
x_j - x_m
end

