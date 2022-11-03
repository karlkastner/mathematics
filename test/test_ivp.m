n=100;
L = 20;
D2 = derivative_matrix_2_1d(n,n,2,'circular');

x = linspace(-L/2,L/2,n)';
y = normpdf(x,0,1);

plot(x,y)

[t,y] = ode23(@(t,y) -D2*y,[0,10],y);

plot(x,y(1:10:end,:)')
legend()
