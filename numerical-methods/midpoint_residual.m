% 2024-09-27 15:22:58.952603374
% Fri 27 Sep 10:50:39 CEST 2024
% Karl Kastner, Berlin
%
%% residual of the implicit midpoint integration rule
%%
%% yc          = (y_next+y)/2)
%% y_next      = y_old + dt*f(yc)
%% res         = y_next - y_old - dt*f(yc)
%% dres/dynext = 1 - dt*df/dyc*dyc/dynext
%% dres/dynext = 1 - 0.5*dt*df/dy_(y=yc)
function [res,Je] = midpoint_residual(fun,dt,x,xold)
	fo = fun(xold);
	if (nargout() < 2)
		[f,Jf] = fun(0.5*x+0.5*xold);
	else
		[f,Jf] = fun(0.5*x+0.5*xold);
		n  = length(x);
		Je = speye(n) - 0.5*dt*Jf;
	end
	res = x - xold - dt*f;
end

