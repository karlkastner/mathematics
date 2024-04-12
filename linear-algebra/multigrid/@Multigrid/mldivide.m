% Mon  8 Jan 13:01:03 CET 2024
function [x,flag,relres,iter,resn] = mldivide(obj,b,x)
	% why 2/3 even though 1 is optimal
	%o = (2/3);
	% number of smoothing steps
	%m = 1;

	% 1 for v-cycle, 2 for w-cycle
	%if (nargin()<9)
	%	gamma = 1;
	%end
	%
	%if (nargin()<7)
	%	maxiter = sum(size(x));
	%end

	%if (nargin()<8)
	%	reltol = sqrt(eps);
	%end
	% implicit euler step
	% dz/dt = a*D2*z + f
	% (I - dt*a*D2) z_k+1 = z_k + dt*f
	
	obj.s(1).x = x;
	obj.s(1).b = b;

	obj.resfun(1);
	resn0 = rms(obj.s(1).res,'all');
	flag  = 0;
	iter  = 0;
	while (1)
		iter = iter+1;
		obj.cycle(1);
		obj.resfun(1);
		resn(iter) = rms(obj.s(1).res,'all');
		if (iter>1 && resn(iter)>=resn(iter-1))
			resn(iter)
			disp('stagnation');
			%break;
		end
		if (resn(iter) <= obj.reltol*resn0)
			break;
		end
		if (iter == obj.maxiter)
			flag = -1;
			disp('maxiter reached');
			break;
		end
	end % while 1
	relres = resn(end)/resn0;
	x = obj.s(1).x;
end % mldivide

