% 2024-09-27 15:22:58.952603374
% Karl Kastner, Berlin
%
%% output
%%	JE : jacobian matrix of the implicit midpoint scheme
%% input 
%% 	jfun : jacobian of the differential equation
%%	dt   : time step
%%	x    : (implicit) value at current time t 
%%	xold : known value of last time step t - dt
function Je = midpoint_implicit_jacobian(jfun,dt,x,xold,returnmat)
	xmid = 0.5*(x+xold);
	Jf = jfun(xmid);
	if (returnmat)
		n  = numel(x);
		Je = speye(n) - 0.5*dt*Jf;
	else
		Je = -0.5*dt*Jf;
		if (size(x,2)==1)
			Je(:,3) = 1 + Je(:,3);
		else
			Je(:,:,5,:) = 1 + Je(:,:,5,:);
		end
	end
end

