% 2016-08-16 21:42:59.640204055 +0200
%
%% least squares by the broyden method
%% for rectangular / non symmetric systems
%%	Numerical  Optimization nocedal
%%	Practical  Methods  of  Optimization fletcher
%% c.f. gerber 1981
%% c.f. fletcher 1978 (more advanced, not used here)
%% c.f. Kelley 1999 ch. 4
%%
%% BGFS:
%% Broyden 1965
%% Fletcher 1970
%% Goldfarb 1970
%% Shanno 1970

function [x nf f H nf_ c x_] = ls_broyden(afun,x,reltol,maxiter,H0,broyden)
	f  = afun(x);
	nf = norm(f)
	nf0 = nf;
	p  = length(x);
	n  = length(f);
	H  = rand(p,n);
	if (~broyden)
		B = H;
	end
	if (~isempty(H0))
		H = H0;
		B = H0;
	end
	% make svd of H small
	%H = 1/(sqrt(p*n))*H;
	H = H;
	nf = inf;
	Hold = H;
	iter = 0;
	nf_(1) = nf0;
	x_ = x;
	while (1)
		fold = f;
		xold = x;
		nfold = nf;
		if (broyden)
			s  = -(H*f);
		else
			s = -(B'\f);
			%s = -( (B'*B) \ B'*f);
		end
		% step, TODO line search
		x  = x+s;
		f  = afun(x);
		nf = norm(f);
%pause
		nf_(iter+2) = nf;
		if (nf < nfold)
			x_ = x;
		end

		if (0) %nf > nfold)
			x = xold;
			f = fold;
			nf = nfold;
			H = Hold;
			break;
		end
		y  = f-fold;
		Hold = H;
		if (broyden)
			% update Hessian
			v  = H'*s/(s'*H*y);
			H = H + (s-H*y)*v';
		try
		c(iter+1) = cond(H);
		catch
		c(iter+1) = NaN;
		end
		else
			% bfgs
		%	H = H + (s'*y + y'*B*y
			h = s;                                                                         
			B = B';
			B = B + (1/(h'*h))*y*h' - (1/(y'*(B*h)))*(B*h)*(y'*B);
			B = B';
		try
		c(iter+1) = 1./cond(B);
		catch
		c(iter+1) = NaN;
		end
		end


		
%		[svd(H)]
%pause
%f
%s
%H*y
%(s'*H*y)
%		nf
%		pause
		iter = iter+1;
		if (iter >= maxiter)
			break;
		end
		if (nf < reltol*nf0)
			break;
		end
	end
end
