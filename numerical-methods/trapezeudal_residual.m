% 2024-09-27 15:22:58.952603374
% Fri 27 Sep 10:50:39 CEST 2024
% Karl Kastner, Berlin
%
%% residual of the (implicit) trapezoidal time step
%%
%% y_next = y_old + dt*f(y_next)
%% res    = y_next - y_old - dt*f(y_next)
%% dres/dynext = 1 - dt*df/dy
function [res,Je] = trapezeudal_residual(fun,dt,x,xold)
	fo = fun(xold);
	if (nargout() < 2)
		[f,Jf] = fun(x);
	else
		[f,Jf] = fun(x);
		n = length(x);
		Je = speye(n) - 0.5*dt*Jf;
	end
	res = x - xold - (0.5*dt)*f - (0.5*dt)*fo;
end

