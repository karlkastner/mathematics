% Mon  4 Sep 18:53:35 CEST 2017
% Karl Kastner, Berlin
%
%% solve intial value problem by the euler forward method
%
% function [y] = ode_euler_forward(odefun,Tspan,y0,opt)
function [t,y] = ode_euler_forward(odefun,Tspan,y0,opt)
	% TODO no magic number
	n = 100;
	dt = (Tspan(2)-Tspan(1))/(n-1);
	t = Tspan(1) + dt*(0:n-1)';

	% dy_dt = zeros(n,1);
	y = zeros(n,1);
	y(1) = y0;
	for idx=2:length(t)
		dy_dt = feval(odefun,t,y(idx-1));
		y(idx) = y(idx-1)+dt*dy_dt;
	end
end

