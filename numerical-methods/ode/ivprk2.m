% Sat 11 Nov 12:01:21 CET 2017
% Karl Kastner, Berlin
%
%% solve initial value problem by the two step runge kutta method
% 
function y = rk2_step(T,x,y,ydot,bcfun)
	% apply bc
	% note: this is not necessary, should be applied in ydot

	% get time step
	dt = T(2)-T(1);

	% TWO step euler vs two step heun
	t = T(1);
	while (t<T(2))
		% limit time step at end
		dt = min(dt,(T(2)-t)*(1+sqrt(eps)));

		% apply bc
		y    = bcfun(t,y);
		%dt_max = 1./max(abs(dydt))

		% step to mid-point
		dydt = ydot(t+dt,x,y_m);
		y_m  = y + 0.5*dt*dydt;

		% apply bc
		y_m  = bcfun(t+0.5*dt,y);

		% step fully
		dydt = ydot(t+0.5*dt,x,y_m);
		y    = y + dt*dydt;
		
		t = t+dt;
	end 
	% check steplength
	%lambda =  dydt./y;
end

