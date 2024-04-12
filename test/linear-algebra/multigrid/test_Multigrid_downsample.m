% 551 2024-01-23 20:58:24.675505318 +0100 test_MG_downsample.m

n = [8,8];
x = rand(n);
a = 1;
d = 1;
L = n;
if (0)
%x = (1:8)'*ones(1,nn(2));
x_ = downsample1(x) 
%x_ = mg.downsample(flat(x),nn);
%reshape(x_,[nn(1)/2,nn(2)])
%reshape(x_,[nn(1),nn(2)/2])'
x_ = mg.downsample1(flat(x'),nn);
x_ = reshape(x_,[nn(1),nn(2)/2])'
end

if (0)
x_ = downsample2(x')' 
%x_ = mg.downsample(flat(x),nn);
x_ = mg.downsample2(flat(x'),nn);
x_ = reshape(x_,[nn(1)/2,nn(2)])'
end

mg_m = Multigrid();
mg_m.init({a},d,L,n,1);

mg_j = javaObject('MG');
mg_j.init(a,d,L,n);


x_m = downsample_2d(x)

x_j = mg_j.downsample12_(x(:)',0);
x_j = reshape(x_j,[n(1)/2,n(2)/2])
x_m - x_j




