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
	iter = 1;
	resold = 0;
	while (1)
pflag = 0;
if (0)
n_ = sqrt(length(x));
h = reshape(x,n_,n_);
clf
imagesc(h)
plot(rms(h-mean(h)))
pause
end
if (pflag)
n_ = sqrt(length(x));
h = reshape(x,n_,n_);
figure(1);
clf
%subplot(2,4,6)
%plot(h0(:,round(end/2),:))
%hold on
plot(h(:,round(end/2)))
%plot(diff(h(:,round(end/2))))
end

		res = fun(x);
		if (pflag)% == iter)
			res_ = reshape(res,n_,n_);
			hold on
			plot(res_(:,round(end/2)))
		%	plot(resold(:,round(end/2)))
drawnow
resold = res_;
		end 
		rmse(iter) = rms(res,'all');
		if (rmse(iter) <= abstol)
			disp('converged');
			break;
		end
		if (iter > maxiter)
			disp('maxiter reached');
			break;
		end
		iter = iter+1;
		switch (method)
		case {'direct'}
			[J] = jfun(x);
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
			
		a = 1;
		% line search
		liter = 0;
		while(1)
			liter = liter+1;
			x_ = x - a*c;

%fun
%dbstop swe_simplified_1d_dzdt
%		[res,bl] = fun(x);
if (0)
%condest(J)
%sum(J(:))
close all
global bl
a
whos x
whos bl
plot([x,res,[diag(J,-1);NaN],diag(J),[diag(J,+1);NaN],-c/30])
legend('x','res','l','c','r')
%x_,-NaN*a*c,res,cvec(bl)])
%NaN*diag(J)])
%imagesc(isfinite(J))
%%plot(isfinite(res))
pause
end
%min(x_(:))
			res_ = fun(x_);
%max(res_)
			if (rms(res_,'all')<rms(res,'all'))
%			x_ = x - 0.5*a*c;
				x = x_;
				break;
			end
			a = a/2;
			if (a < sqrt(eps))
				disp('not a descend direction')
				iter = maxiter+1;
				break; %error('here')
			end
		end
if(pflag)
rms(res_,'all')
	c_ = reshape(c,n_,n_);
			hold on
			plot(-c_(:,round(end/2)))
			plot(-a*c_(:,round(end/2)))
pause()
end
a

	end

function x = mgmfun(rhs)
	rhs = reshape(rhs,n(1),n(2),nvar);
	x = zeros(n(1),n(2),nvar);
	%reshape(x,n(1),n(2),nvar);
	x = mg.cycle1(rhs,x);
	x = flat(x);
end

end

