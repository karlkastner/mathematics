% Mon  8 Jan 18:12:06 CET 2024
%
% c.f. hackbusch 
function [x,resn]  = mg_heat_2d_simple(a0,a,b,x,n,dt,dx,reltol,maxiter)
	o = 2/3;
	m = 1;
	%D2 = derivative_matrix_2_1d(n,(n-1)*dx,2,'circular');
	L = (n-1)*dx;
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,{'circular','circular'});
	D2 = D2x+D2y;
	%D2 = derivative_matrix_2_1d(n,(n-1)*dx,2,'hdirichlet');
	I  = speye(prod(n));
	A  = a0*I - (a*dt)*D2;

	A_C    = {A};
	iD_C   = {};
	down_C = {};
	dt_C  = {};
	up_C   = {};
	ut_C = {};
	k  = 1;
	n_ = n(1); %length(b);
	while (1)
		iD_C{k}   = 1./diag(A_C{k});
		if (n_>=2)
			[down_C{k},dt_C{k}] = downsampling_matrix_2d([n_,n_]); %,'mg');
			[up_C{k},dt_C{k}]   = upsampling_matrix_2d([n_,n_]);
			I_ = speye(n_);
			%A_C{k+1}  = kron(down_C{k},dt_C{k})*A_C{k}*kron(up_C{k},ut_C{k});
test with dsm/usm
????
			A_C{k+1}  = kron(down_C{k},I_)*kron(I_,dt_C{k})*A_C{k}*kron(I_,dt_C{k})*kron(down_C{k},I_);
			n_        = n_/2;
			k         = k+1;
		else
			break;
		end
	end
	iter   = 0;
        x = x(:);
	b = b(:);
	resn0  = rms(A*x - b,'all');
	resn   = [];
	dresn_ = [];
	while (1)
		iter       = iter+1;
		x          = v_cycle(1,b,x);
		res        = resfun(1,b,x);
		resn(iter) = rms(res,'all');
		if (resn(iter) <= reltol*resn0)
			break;
		end
		if (iter == maxiter)
			warning('maxiter reached');
			break;
		end
	end
	x = reshape(x,n);

	function res = resfun(k,b,x)
		res = (A_C{k}*x - b);
	end

	function [x] = v_cycle(k,b,x)
	if (size(x,1) >= 2)
		% pre-smooth, hackbusch 2.5.2.b
		for idx=1:m
			% since x is zero except in the first level, the first mvec can be saved
			res = resfun(k,b,x);
			x = x - o*(iD_C{k}.*res);
		end

		% recurse
		% note that it is necessary to restrict the defect, not the solution,
		% since only smooth functions can be well represented on coarser grids

		% compute residual 2.5.2c
		res  = resfun(k,b,x);

		% downsample (restrict) 2.5.2c
		n = sqrt(numel(x));
		res_ = flat(down_C{k}*reshape(res,n,n)*dt_C{k});
		x0   = zeros(size(res_));
		% recurse to approximated error
		e_   = v_cycle(k+1,res_,x0);
		% 2.5.2c4
		e = flat(up_C{k}*reshape(e_,n/2,n/2)*ut_C{k});
		%e    = up_C{k}*e_;
		% correct
		x    = x - e;

		% post-smooth
		for idx=1:m
			res = resfun(k,b,x);
			x   = x - o*(iD_C{k}.*res);
		end

	else
%		x   = x + o*(iD_C{k}.*res);
		x = A_C{k} \ b(:);
		%x = reshape(x,size(b));
	end
end

end % mg_heat_1d_simple

