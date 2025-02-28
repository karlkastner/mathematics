% 2024-11-08 20:29:42.329889707 +0100
% Karl Kastner, Berlin
function [x,rmsr,C,r2,search_direction] = gauss_newton_(resfun,x,tol,maxiter)
	if (nargin()<3 || isempty(tol))
		tol     = 1e-3;
	end
	if (nargin()<4)
		maxiter = 10;
	end
	x       = cvec(x);
	nres    = inf;
	iter    = 0;
	while (1)
		iter = iter+1;
		[res,J] = resfun(x);
		rmsr = rms(res);
		
		grad = J'*res;
		JtJ = J'*J;
		search_direction = -JtJ\grad;
		% line search
		% TODO this does not check if the gradient decreases
		a = 1;
		liter = 0;
		while (1)
			x_ = x - a*search_direction;
			res_ = resfun(x_);
			rmsr_ = rms(res_);
			if (rmsr_ < rmsr)
				rmsr = rmsr_;
				x = x_;
				break;
			end
			a = a/2;
			if (a<sqrt(eps))
				disp('not a descend direction');
				iter = maxiter+1;
				break;
			end
		end
		if (grad'*grad < tol^2)
			break;
		end
		if (iter>maxiter)
			disp('maxiter reached');
			break;
		end
	end
	C = rmsr^2*inv(JtJ);
end

