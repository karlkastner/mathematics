% 2021-06-23 16:03:48.389683349 +0200
function [t,y] = oscillator_noisy(T,y0,c,s)
	%ydot = @(t,x) [0,-1; 1, 0]*x + (1-hypot(x(1),x(2))).*x
	%ydot = @(t,x) [0,-1; 1, 0]*x + (1-hypot(x(1),x(2))).*x
	%ydot = @(t,x) [0,-1; 1, 0]*x + c*(1-hypot(x(1),x(2)).^2).*x + s*randn(2,1);
	ydot = @(t,x) ( [0,-1; 1, 0]*x ... % oscillator
                        + s*randn(2,1) ... % perturbation
        		+ c*(1-hypot(x(1),x(2)).^2).*x % amplitude recovery (non-linear)
			);
	% c*(x - x^3 - y^2*x)
        % + c*(x-x.^2) <- simpler, but unstable
	if (0 == s)
		[t,y] = ode23s(ydot,T,y0);
	else
		%[t,y] = euler(ydot,T,y0);
		dt = 1e-4;
		[t,y] = ivp_euler_forward2(ydot,T,y0,dt);
	end
end

