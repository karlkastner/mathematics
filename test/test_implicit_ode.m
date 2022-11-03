% note : the solver in mass form does not take fewer steps
% in this case, the jacobian with mass matrix is simpler
function [t,y]=test_implicit()

	% -y/(1+y)
	y0 = 1;
	T = [0,10];
	p = 2;

	solver = @ode23t
%	solver = @ode15s
%	solver = @ode23tb
	
	opt = odeset('Jacobian',@(t,y) -2*y/(1+y) + y^2/(1 + y)^2);
	[t,y] = solver(@(t,y) - y^2/(1 + y),T,y0);
	length(t)
	
	clf
	plot(t,y,'.')
	opt = odeset('Mass',@(t,y) (1+y),'Jacobian',@(t,y)-2*y);
	[t,y] = solver(@(t,y) - y^2,T,y0,opt);
	hold on
	plot(t,y,'o')
	length(t)

        opt = odeset('Jacobian',@(t,y) -2*y/(1+y) + y^2/(1 + y)^2);             
	[t,y] = trapezoidal_fixed(@(t,y) - y^2/(1 + y), t, y0,opt)
	plot(t,y,'-')
end

