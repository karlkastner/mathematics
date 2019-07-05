% Thu 16 Nov 13:09:39 CET 2017
% Karl Kastner, Berlin
%
%% single trapezoidal step
function y = step_trapezoidal(t,y0,ydot,dt)
	abstol = sqrt(eps);

	% derivative at origin
	dy0  = ydot(t,y0);
	dy   = dy0;
	yold = y0;

	% single step
	while (true)
		y = y0 + 1/2*dt*(dy0 + dy);
		if (max(abs(y-yold)) < abstol)
			break;
		end

		dy   = ydot(t,y);
		yold = y;
	end % while true
end

