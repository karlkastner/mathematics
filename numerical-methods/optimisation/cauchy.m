% 2014-08-25 19:03:25.598490341 +0200
function	[x flag F] = cauchy(func,x)
	tol = sqrt(eps);
	maxj = 12;
	imax = 10000;
	l = 0.5;
	idx = 0;
	fold = 0;
	while (1)
		idx = idx+1;
		if (idx > imax)
			'imax'
			flag = 0;
			return;
		end
		f0 = feval(func,x);
		F(idx,1) = f0;
		g = grad(func,x);		
		%if (norm(f0) < tol*norm(x))
		G(idx,1) = norm(g);
		if (abs(fold - f0) < tol)
			fold-f0
		%if (norm(g) < tol*norm(x))
			'converged'
			plot([F G])
			pause
			flag = 1;
			break;
		end
		l = 2*l;
		jdx = 0;
		while (feval(func,x-g*l) >= f0)
			l = 0.5*l;
%			[f0 feval(func,x-g*l)]
%			[ l idx]
%			plot(F)
			jdx = jdx+1;
			if (jdx > maxj)
				'maxj'
				flag = 0;
				return;
			end
		end
		x = x-g*l;
		fold = f0;
	end
end

