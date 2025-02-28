% 2024-09-27 10:40:20.096248924 +0200
% Karl Kastner, Berlin
%
%% Gauss-Newtown solver, a.o. for solving implicit time_steps of nonlinear PDEs
function [x,rmse,res,iter,out] = gauss_newton(fun,x,tol,maxiter,method,jfun,L,n,nvar,initarg,jfun2)
	if (nargin()<3 || isempty(tol))
		 abstol = 1e-7;
		 reltol = abstol;
	end
	if (nargin()<4)
		maxiter = 10;
	end
	x       = cvec(x);
	iter = 0;
	resold = 0;
	while (1)
		iter = iter+1;

		res = fun(x);
		rmse(iter) = rms(res,'all');
		if (rmse(iter) <= abstol)
			disp('converged');
			break;
		end
		if (iter > maxiter)
			disp('maxiter reached');
			break;
		end
		switch (method)
		case {'direct'}
			J = jfun(x);
			c =  (J \ res);
			out.J = J;
		case {'ilu-ks'}
			J = jfun(x);
			[L,U] = ilu(J);
			n = length(J);
			%c = bicgstabl(J,res,reltol,n,L,U);
			c = gmres(J,res,10,reltol,n,L,U);
			out.J = J;
		case {'mg'}
			mg = Multigrid();
			xi = x;
			xi(:,:,end+1:end+size(initarg,3)) = initarg;
			mg.init_fun(jfun,L,n,nvar,xi);
			mg.opt.reltol = 0.1*reltol;
			mg.opt.abstol = 0.1*reltol;
			x0 = zeros(n(1),n(2),nvar);
			[c,flag,relres,mgiter,resn] = mg.mldivide(reshape(res,n),x0);
			%c = c(:);
			out.mg = mg;
			% mg.check()
		case {'mg-ks'}
			mg = Multigrid();
			xi = reshape(x,n);
			xi(:,:,end+1:end+size(initarg,3)) = initarg;
			mg.init_fun(jfun2,L,n,nvar,xi);
			J = jfun(x);
			x_ = mgmfun(res);
                        c = bicgstabl(J,res,reltol,[],@mgmfun);
		otherwise
			disp(method);
			error('unknown method');
		end
			
		% line search
		a = 1;
		liter = 0;
		while(1)
			liter = liter+1;
			x_ = x - a*c;
			res_ = fun(x_);
			rms_ = rms(res_,'all');
			if (rms_<rmse(iter)) % (res,'all'))
				rmse(iter) = rmsr_;
				x = x_;
				break;
			end
			a = a/2;
			if (a < sqrt(eps))
				disp('not a descend direction');
				iter = maxiter+1;
				break;
			end
		end

	end

function x = mgmfun(rhs)
	rhs = reshape(rhs,n(1),n(2),nvar);
	x = zeros(n(1),n(2),nvar);
	%reshape(x,n(1),n(2),nvar);
	x = mg.cycle1(rhs,x);
	x = flat(x);
end

end

