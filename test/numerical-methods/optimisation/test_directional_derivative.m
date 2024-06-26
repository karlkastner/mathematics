% 2016-08-22 14:52:05.901387964 +0200
fun   = @(x) exp(-2*x)
order = 4;
x0    = 0;
dir   = 1;
h     = eps^0.125;
directional_derivative(fun,x0,dir,h,order)

