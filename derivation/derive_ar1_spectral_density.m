% 2021-06-26 20:16:32.236412500 +0200

syms x t k r
%assume(k > 0)
assume(in(t,'real'))

% from acf to ft  : ft = sum_k exp(2i pi k x) a dk
% from ft  to acf : a = int exp(2i pi k x) ft dx

if (0)

% acf for lp:
a_lp = exp(-k*abs(x))

% acf for all-pass
a_ap = 1;

% 1(x-t)*exp(-k*x)
int(exp(-k*x),x,0,inf)
int(exp(-k*(t-x)),x, -inf,0)
int(exp(-k*(t+x)),x, 0,inf)
end

%a = r^-abs(k)
%int(exp(2*i*pi*x*k)*a,x)
symsum(exp(2*i*pi*x*k)*r^-abs(k),k)
% piecewise(r*exp(pi*x*2i) == 1, 0, r*exp(pi*x*2i) ~= 1, 1/(r*exp(pi*x*2i) - 1))

% wolfram: sum exp(2*i*pi*x*k)*r^-abs(k), k = -infty,  infty
% f = (1 - r^2)/((1 - r*exp(-2i*pi*x))*(1 - r*exp(2*i*pi*x)))
%f = (1 - r^2)/(1 - r*(exp(-2i*pi*x)+exp(2i*pi*x)) + r^2)
f2 = (1 - r^2)/(1 - 2*r*cos(2*pi*x) + r^2)
% inverse for test : int exp(2*i*pi*x*k)*(1 - r^2)/(1 - r*(exp(-2i*pi*x)+exp(2i*pi*x)) + r^2) dx, x = -infty, infty
% (computation time exceeded)

% acf of high-pass
% wolfram: sum exp(-2*i*pi*x*k)*(1-(1 - r^2)/(1 - 2*r*cos(2*pi*x) + r^2)), k = -infty, infty
% wolfram: sum exp(-2*i*pi*x*k)*(1 - r^2)/(1 - r*(exp(-2i*pi*x)+exp(2i*pi*x)) + r^2), k = -infty, infty


% sum exp(-2*i*pi*x*k)*((1 - r^2)/(1 - r*(exp(-2i*pi*x)+exp(2i*pi*x)) + r^2)*(1-(1 - r^2)/(1 - r*(exp(-2i*pi*x)+exp(2i*pi*x)) + r^2))), k = -infty, infty
