nx = 100;
L = 10;
w = 1;
x = linspace(-L/2,L/2)';

a = 1;
t = 0.09;
y0 = abs(x)<w;

a = -1;

k = 20;
y = y0;
for idx=1:k
	y = advect_analytic(t,y,L,a);
end

D1 = derivative_matrix_1_1d(nx,L,-sign(a),'circular');

ode = @(t,y) a*D1*y;
[t,yy] = ode23(ode,[0,k*t],y0);
y(:,2) = yy(end,:);

plot(x,[y0,y]);
