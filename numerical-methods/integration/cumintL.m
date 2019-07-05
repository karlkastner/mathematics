% Wed  1 Nov 11:34:07 CET 2017
%% cumulative integral from left to right
% Note: a previous version used to give x as the second argument
% input:
% y  : [n x m]
% dx : scalar or [n-1 x m]
% 
% output : int_dy_dx : [nxm]
function inty = cumintL(y,arg2);
	transposed = false;
	if (isvector(y))
		if (isrow(y))
			transposed = true;
		end
		y  = cvec(y);
	end
	if (isscalar(arg2))
		dx = arg2;
	else
		if (isvector(arg2))
			arg2 = cvec(arg2);
		end
		dx = diff(arg2);
	end
	nx = size(dx,1);
	ny = size(y,1);
	if (isscalar(dx) || nx == ny-1)
	% values given at end-points
	% integration by the trapezoidal rule
	inty = [ zeros(1,size(y,2)); ...
                 0.5*cumsum(bsxfun(@times,dx,(y(1:end-1,:)+y(2:end,:)))) ];
	else
	% values given at mid-points
	% integration by the midpoint rule
	inty = [ zeros(1,size(y,2));
		 cumsum(bsxfun(@times,dx,y)) ];
	end

	if (transposed)
		inty = inty.';
	end
end

