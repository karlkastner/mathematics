% 2016-08-15 17:09:00.677048425 +0200
if (0)
% sample, not parameter
x  = rand(100,1);
c0 = randn(5,1)
%c0(1) = 0;
%c0(3) = 0;
%c0(end) = 100;
%f = (x-x0).^4;
fun  = @(c) polyval(c,x) - polyval(c0,x)
%polyval(x,c)

c0_ = 2*ones(size(c0));
dc  = 0.1;
[c_ f_] = dud(fun,c0_,dc)
[c0 c_]
end

%x  = [-1.2,1]
%dx = 0.1*x;
%dud(@rosenbrock,x,dx)

x = (1:10)'/10;

%c = [1 1];
%c = [5;20];
%c = [2.5,10];
c = [1;2;3];
dc = [0.1,0.2,-0.1]; %[1;2];
dud(@(c) box2(c,x),c,dc)
