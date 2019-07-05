% 2015-06-22 08:41:29.485788940 +0200
% Karl Kastner, Berlin
%% fit a polynomial function
function [param A] = fit_(X,Y,W,order)
	A      = vander_1d(X,order);
	if (nargin() > 3 && ~isempty(W))
		Y = double(Y);
		A = double(A);
		if (1 == size(W,3))
			if (isvector(W))
				W = diag(sparse(W));
			end
			param = (A'*W*A) \ (A'*W*Y);
		else
			for idx=1:size(W,3)
				wi = W(:,:,idx);
				if (isvector(wi))
					wi = diag(sparse(wi));
				end
				param(:,idx) = (A'*wi*A) \ (A'*wi*Y(:,idx));
			end
		end
	else
		% solve by QR factorisation
		[Q R]  = qr(A,0);
		%param  = R \ (Q'*Y);
		opts.UT = true;
		param = linsolve(R,Q'*Y,opts);
	end
end % fit_

