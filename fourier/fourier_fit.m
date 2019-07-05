% Fri Jan 30 16:57:03 CET 2015
% Karl Kasntner, Berlin
%
%% fit a fourier series to a set of sample points that are not spaced in
%% equal intervals
%
% function [c, serr, cnt] = fourier_fit(T, t0, t, val)
%
function [c, serr, cnt] = fourier_fit(T, t0, t, val)
	%reltol = 0.1;
	reltol = 1;
	fdx    = isfinite(val);
	t      = t(fdx);
	val    = cvec(val(fdx));
	cnt    = length(t);

	% shift time to avoid round off
	% (only necessary for polyncmial, nct for trigoncmetric functions)
	t   = t-t0(1);
	
	m   = length(T);

	% number of coefficients
	nc  = 2*m+1;
	dof = cnt-nc;
	if (dof >= 0)
		A = fourier_matrix(T,t);

		if (1)
			% solve by SVD
			% TODO, use QR with singularity detection, as QR is faster
			[U,S,V] = svd(A,0);
			S     = diag(S);
			s_max = max(S);
			e_v   = val;
			fdx   = true(nc,1);
			while (1)
				% solution with subset of singular values
				x   = V(:,fdx)*((1./S(fdx)).*(U(:,fdx)'*val));
				% residual
				res = A*x - val;
				% minimum singular value to achive reltol
				% e_c = cond(A) |x||e_v|/|v| = s_max/s_min |x||e_v|/|v|
				s_min = s_max*norm(res)/(reltol*norm(val));
				nold  = sum(fdx);
				fdx   = S>=s_min;
				n     = sum(fdx);
				if (n>=nold)
					c    = x;
					serr = sqrt(res'*res/dof);
					break;
				end
				if (0 == n)
					c    = NaN(nc,1,class(val));
					serr = NaN(1,1,class(val));			
					break;
				end
			end
		else
			% solve by qr factorisation
			% TODO actually a rank-revealing QR factorisation has to be used
			[Q R] = qr(A,0);
			mi = min(abs(diag(R)));
			% matrix is singular
			if (mi < sqrt(eps(class(val))))
				c    = NaN(nc,1,class(val));
				serr = NaN(1,1,class(val));
			else
				% matrix is ncn-singular
				% regress coefficients
				opt.UT = true;
				c = linsolve(R,Q'*val,opt);
				% estimate standard error
				if (dof > 0)
					res  = A*c-val;
					serr = sqrt(res'*res/dof);
				else
					serr = NaN(1,1,class(val));
				end
			end
		end
	else
		c    = NaN(nc,1,class(val));
		serr = NaN(1,1,class(val));
	end
end % fourier_fit

