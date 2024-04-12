% Mon  8 Jan 13:01:03 CET 2024
function [x,resn] = mg_heat_step_1d(a0,a,b,x,dt,dx)
	o = 2/3;
	m = 1;
	% implicit euler step
	% dz/dt = a*D2*z + f
	% (I - dt*a*D2) z_k+1 = z_k + dt*f
	reltol = 1e-7;
	maxiter = 100;

	% TODO only for x0 = 0
	resn0 = rms(b);
	iter = 0;
	while (1)
		iter    = iter+1;
		x   = v_cycle(b,x,dx);
		res = resfun(b,x,dx);
		resn(iter) = rms(res);
		if (resn(iter) <= reltol*resn0)
			break;
		end
		if (iter == maxiter)
			warning('maxiter reached');
			break;
		end
	end

	function res = resfun(b,x,dx)
		% R = -dt*a/dx^2*[1,0,1]
		%  r = b - (a0*x - dt*a*[1,-2,1]/dx^2*x)
		res = ( (  ( (a0+2*dt*a/dx^2))*x ...
			    - (dt*a/dx^2)* ...
			      (   circshift(x,-1) ... 
			        + circshift(x,+1) ...
			      ) ...
			   ) - b ...
		     );
	end % resfun

	function x = jacobi_step(b,x0,dx_)
		% (D+R)*x = b
		% x =     o*D^-1*(b-R*x);
		% x = x + o*D^-1*(b-A*x);

		d = a0 +  dt*a*2/dx_^2;

		res = resfun(b,x0,dx_);

		x = x0 - o*(res./d);
	end

	function [x] = v_cycle(b,x,dx)
		if (length(x) > 1)
			% pre-smooth, hackbusch 2.5.2.b
			for idx=1:m
				x = jacobi_step(b,x,dx);
			end	

			% n.b. for the w-cycle, this is simply repated twice
			% residual
			res = resfun(b,x,dx);

			% downsample (restrict)
			res_ = downsample_1d(res);

			% recurse to approximate error
			x0 = zeros(length(b)/2,1);
			e_  = v_cycle(res_,x0,2*dx);

			% upsample and correct
			x   = x - upsample_1d(e_);

			% post-smooth
			for idx=1:2
				x = jacobi_step(b,x,dx);
			end
		else
			% (2*(a0-dt*a/dx^2))*x = b
			%x = b./(a0+dt*a*2/dx^2);
			x = b./(a0); %+dt*a*2/dx^2);
		end
	end % v-cycle		
end

