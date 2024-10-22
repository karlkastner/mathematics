% 2024-09-27 15:22:58.952603374
% Karl Kastner, Berlin
%% jacobian of the implicit euler time step
%% requires function determining the jacobian of the differential equation
function Je = euler_implicit_jacobian(jfun,dt,x,returnmat)
	Jf = jfun(x);
	if (returnmat)
		n  = numel(x);
		Je = speye(n) - dt*Jf;
	else
		Je = -dt*Jf;
		if (size(x,2)==1)
			Je(:,3) = 1 + Je(:,3);
		else
			Je(:,:,5,:) = 1 + Je(:,:,5,:);
		end
	end
end

