function [x,flag,resn,resn0] = mldivide_jacobi(A,b,x,reltol,maxiter)
	o = 2/3;
	flag = 0;
	if (nargin()<4)
		reltol = sqrt(eps);
	end
	if (nargin()<5)
	maxiter=1e3;
	end
	tol = 1e-4;
	res = A*x-b;
	resn0 = rms(res);
	D = full(diag(A));
	iter=0;
	while (1)
		iter = iter+1;
		% b=A*x
		% b = D*x + R*x
		%x = (b-R*x)./D;
		%x = (b-R*x-Dx+Dx)./D;
		x = x - o*(A*x-b)./D;
		res = A*x-b;
		resn = rms(res);
		if (resn <= tol*resn0)
			break
		end
		if (iter > maxiter)
			flag = -1;
			break;
		end
	end
end
