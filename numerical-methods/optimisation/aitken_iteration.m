% Sun 13 Jun 21:51:31 CEST 2021
% aitken fixed point interation
%
% Note : this should work for flatted matrices just fine, as the inner products
%        are transpose-invariant
% TODO for matrix : messaoudi 1997, matrix extrapolation                                            
%              messaoudi 1993                                                                  
%              jilou 2015 
%
% function [x0, cflag, iter] = aitken_iteration(fun, x0, opt,flag)
function [x0, cflag, iter] = aitken_iteration(fun, x0, opt,flag)
	if (nargin()<3)
		opt = struct();
	end
	if (~isfield(opt,'abstol'))
		opt.abstol  = 1e-7;
	end
	if (~isfield(opt,'reltol'))
		opt.reltol = 1e-3;
	end
	if (~isfield(opt,'maxiter'))
		opt.maxiter = length(x0);
	end
	if (~isfield(opt,'relaxation'))
		opt.relaxation = 1;
	end
	if (~isfield(opt,'verbose'))
		opt.verbose = 0;
	end
	if (~isfield(opt,'miniter'))
		opt.miniter = 0;
	end
	if (nargin()<4)
		flag = true;
	end
	iter = 0;
	cflag = 0;
	while (1)
		iter = iter+1;
		x1   = fun(x0);
		x2   = fun(x1);
		D2x  = (x2-x1)-(x1-x0);
		if (~flag)
			D1x  = (x2-x1);
			dx = D1x.^2./D2x;
			x3 = x2 - dx; 
		else
			% c.f. chow 1984, citing jennings

			% for vectors
			% sidi 2017
			% dx = d2x \ (d1x.^2)
			D1x = x1 - x0;
			D1x21 = x2-x1;                                 
			% D2x = (x2 - x1) - (x1 - x0);                                                   
		%	dx = (D1x.'*D1x)./(D1x.'*D2x)*D1x;  % sidi 2017
			lambda = (D1x21.'*D1x)./(D1x21.'*D2x)
			dx = lambda*D1x; % jennings
			% note that indeed x0 - dx converges faster than x1-dx or x2-dx
			x3 = x0 - dx; 
		end

		x0   = x3;
		ndx  = norm(dx);

		% check for convergence
		if (ndx < opt.reltol*norm(x0))
			cflag = 1;
			break;
		end

		%if (iter > 1 && ndx < opt.reltol*ndxold)
		%	cflag = 2;
		%	break;
		%end
		if (iter > opt.maxiter)
			warning('stopped before converging because maximum number of iterations was reached');
			break;
		end
	end % while 1
end % aitken_iteration

