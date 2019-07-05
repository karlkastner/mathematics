% Sun 12 Jun 16:48:27 CEST 2016
% Karl Kastner, Berlin
%% optimize by quadratic programming
function [x f g] = quadratic_programming(fun,x,opt)
	if (nargin() < 3)
		opt = struct();
	end
	if (~isfield(opt,'MaxIter'))
		opt.MaxIter = 20;
	end
	if (~isfield(opt,'abstol'))
		opt.abstol = 1e-7;
	end	
	if (isfield(opt,'lambda'))
		lambda = opt.lambda;
	else
		% initial levenberg-marquadrd parameter
		lambda = 10;
	end

	maxlm = 100;
	idx = 0;
	f   = inf;
	fopt = inf;
	while (true)
		idx=idx+1;
		fold = f;
		% stop if maximum number of iterations has been reached
		if (idx > opt.MaxIter)
			printf('Stopped at iteration %d without convergence as maximum number of iterations has been reached\n',idx);
			break;
		end
		% compute value and derivatives of local optimisation function
		% not necessary to reduce points to local coordinates here
		[f g H] = fun(x);
%		[f g H] = obj.objective_angle(Pl(:),Tl,[]);
%		[H_ g_ f_] = hessian(@(Pl) obj.objective_angle(reshape(Pl,[],2),Tl), Pl(:));
		if (f >= fold - opt.abstol)
			printf('Stopped at iteration %d as change in function value fell below absolute tolerance\n',idx);
			break;
		end

		% levenberg-marquardt iteration
		jdx=0;
		while (1)
		jdx = jdx+1;

		% compute linearized optimum
		if (~issparse(H))
			dx = (H+lambda*eye(size(H)))\g;
		else
			dx = minres(H+lambda*speye(size(H)),g,[],length(g));
		end
		xopt = x - dx;
%		idx
%		lambda
%		max(abs(dx))
		f
		fold
		fopt

		% new function value
		fopt = fun(xopt);
		fopt
		sum(xopt)

		if (isfinite(fopt) && fopt < f)
			lambda = 0.5*lambda;
			break;
		end
		if (jdx > maxlm)
			fprintf('No improvement during lm iteration');
			return;
		end

		% increase lm-regularisation
		
%pause
		lambda = 2*lambda

		end % while

		% perform line search
		%h0 = 1;
		%[fopt Popt_ h] = line_search(fun, xopt-x, 0, h0, [], [], maxit);
		%if (0)
		%	% in rare cases a point can move further than it's closest neighbour and the triangulation may become inconsistent
		%	% therefore the step is limited to halve the distance to the closesest neihbour
		%	dx = bsxfun(@minus,Pl(fdx,1),Pl(fdx,1)');
		%	dy = bsxfun(@minus,Pl(fdx,2),Pl(fdx,2)');
		%	l  = sqrt(dx.^2+dy.^2) + 1e7*eye(size(dx));
		%	lmin = min(l)';
		%	delta = hypot(Popt(:,1)-Pl(fdx,1),Popt(:,2)-Pl(fdx,2));
		%	h = min(0.5*min(lmin./delta),h);
		%end
		% step
		%xopt(fdx) = (1-h)*x(fdx) + h*xopt(fdx);
		x    = xopt;
	end % while true
end % quadratic_programming

