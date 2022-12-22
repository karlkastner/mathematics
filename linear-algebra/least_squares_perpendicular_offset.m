% Tue 11 Oct 09:45:13 CEST 2022
% least squares line fit minimizing squared deviation to lines,
% this assumes equal uncertainty in x and y
% compated to uncertainty only in y for ordinary least squares
% c.f. Sardelis, D. and Valahas, T., Least Squares Fitting--Perpendicular Offsets, Wolfram
function [slope, intercept,res2,yp] = least_squares_perpendicular_offset(x,y,w)
	x = cvec(x);
	y = cvec(y);
	if (nargin()<3)
		w = ones(size(x,1),1);
	end
	w = w/sum(w);
	% B = 0.5*( mean(y.^2) - mean(y)^2 - mean(x.^2) + mean(x).^2)./(mean(x).*mean(y) - mean(x.*y));
	B  = 0.5*( (w'*y.^2) - (w'*y)^2 - (w'*x.^2) + (w'*x).^2)./((w'*x).*(w'*y) - w'*(x.*y));
	slope = -B + sqrt(B^2 + 1);
	slope(2) = -B - sqrt(B^2 + 1);
	if (0 == (w'*x).*(w'*y) - w'*(x.*y))
		slope = [0,0];
	end

	intercept = (w'*y) - slope*(w'*x);
	yp = [intercept+x*slope];
	% squared residual
	res2 = (w'*(y - yp).^2)./(1+slope(1).^2);

	% the root has two solutions but only one is optimal
	% select the optimal solution
	if (res2(1)<res2(2))
		intercept = intercept(1);
		slope  = slope(1);
		res2 = res2(1);
		yp = yp(:,1);
	else
		intercept = intercept(2);
		slope = slope(2);
		res2 = res2(2);
		yp = yp(:,2);
	end
end % ls_perpendicular_offset

