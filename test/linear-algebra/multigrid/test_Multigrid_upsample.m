% 2024-01-23 20:57:39.382772963 +0100

n = [8,8];
x = rand(n/2);
a = 1;
d = 1;
L = n;


if (0)
x_m = upsample1(x)
%x_j = javaMethod('MG','upsample1',flat(x'),n);
x_j = mg.upsample2(flat(x),n);
x_j = reshape(x_j,[2*n(1),n(2)])
x_j - x_m
end

if (0)

x_m = upsample2(x)
x_j = mg.upsample1(flat(x),n);
%x_j_ = reshape(x_j,[2*n(1),n(2)])
x_j = reshape(x_j,[n(1),2*n(2)])
x_j - x_m
end

mg_m = Multigrid();
mg_m.init({a},d,L,n,1);

mg_j = javaObject('MG');
mg_j.init(a,d,L,n);
x_m = upsample_2d(x)

x_j = mg_j.upsample12_(flat(x)',0);
%x_j_ = reshape(x_j,[2*n(1),n(2)])
x_j = reshape(x_j,[n(1),n(2)])
x_j - x_m

