% Mon  8 Jan 13:01:03 CET 2024
%% solve linear system using multigrid
function [x,flag,relres,iter,resn] = mldivide(obj,b,x)

	obj.s(1).x = x;
	obj.s(1).b = b;

	obj.resfun(1);
	rmsb  = rms(b,'all');
	%resn0 = rms(obj.s(1).res,'all');
	%scale = rms(x,'all');
	flag  = 0;
	iter  = 0;
	while (1)
		iter = iter+1;
		obj.cycle(1);
		obj.resfun(1);
		resn(iter) = rms(obj.s(1).res,'all');
		if (iter>1 && resn(iter)>=resn(iter-1))
			disp('stagnation');
			break;
		end
		if (   resn(iter) <= obj.opt.reltol*rmsb ...
		    || resn(iter) < obj.opt.abstol ...
		   )
			break;
		end
		if (iter == obj.opt.maxiter)
			flag = -1;
			disp('maxiter reached');
			break;
		end
	end % while 1
	relres = resn(end)/rmsb; %resn0;
	x = obj.s(1).x;
end % mldivide

