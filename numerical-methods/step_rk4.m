% Tue 12 Mar 15:07:50 CET 2024
function y = step_rk4(t0,dt,y0,ydot)
	a = dt*ydot(t0, y0);
	b = dt*ydot(t0+dt/2,y0 + a/2);
	c = dt*ydot(t0+dt/2,y0 + b/2);
	d = dt*ydot(t0+dt,y0 + c);
	y = y0 + 1/6*(a + 2*b + 2*c + d);
end


