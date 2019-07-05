% Fri Jan 30 16:57:03 CET 2015
% Karl Kasntner, Berlin
%
%% fit a continous fourier series to a set of sample points not sampled
%% at equal intervals
%
function [c, serr, cnt] = fourier_regress(T, t0, t, val)
% TODO maximum frequency can be determined by maximum gap length
% max(diff([t(end)-T,t]);

	fdx = isfinite(val);
	t   = t(fdx);
	% shift t to avoid round off (only necessary for polynomial)
	t   = t-t0(1);
	val = val(fdx);
	cnt = length(t);
	
	m   = length(T);
	no  = 2*m+1;
	dof = cnt-no;
	if (dof >= 0)
		% usually the first fourier coefficient has the pre-factor 1/2
		% here it is left with 1 so that it equals the mean (in case frequencies are orthogonal)
		A = zeros(cnt,no);
		A(:,1) = 1;
		for idx=1:m
			A(:,2*idx)   = sin(2*pi*t/T(idx));
			A(:,2*idx+1) = cos(2*pi*t/T(idx));
		end
		% solve by qr factorisation
		[Q, R] = qr(A,0);
		% TODO use rank revealing QR factorization
		mi = min(abs(diag(R)));
		% matrix is singular
		if (mi < 1e-7)
			c    = NaN(no,1,class(val));
			serr = NaN(1,1,class(val));
			return;
		end
		% regress coefficients
		%c   = A \ val;
		opt.UT = true;
		c = linsolve(R,Q'*val,opt);
		% estimate standard error
		if (dof > 0)
			res  = A*c-val;
			serr = sqrt(res'*res/dof);
		else
			serr = NaN(1,1,class(val));
		end
	else
		c    = NaN(no,1,class(val));
		serr = NaN(1,1,class(val));
	end
end % fourier_regress

