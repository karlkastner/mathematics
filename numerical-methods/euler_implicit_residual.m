% 2024-09-27 10:50:39 CEST
% Karl Kastner, Berlin
%
%% residual of the implicit-euler time step
%%
%% y_next = y_old + dt*f(y_next)
%% res    = y_next - y_old - dt*f(y_next)
%% dres/dynext = 1 - dt*df/dy
function [res,Je] = euler_implicit_residual(fun,dt,x,xold,nvar)
	if (nargout() < 2)
		[f,Jf] = fun(x);
	else
		[f,Jf] = fun(x);
		n = length(x);
		Je = speye(n) - dt*Jf;
	end
	res = x - xold - dt*f;
end

