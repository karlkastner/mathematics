% 2018-08-29 14:27:21.089395412 +0200
%% integrate y along x with the trapezoidal rule
function cumint_y = int_trapezoidal(x,y)
	if (isvector(y))
		y = cvec(y);
	end
	if (isvector(x))
		x = cvec(x);
	end
	dx = diff(x);
	cumint_y = [zeros(1,size(y,2)); 0.5*cumsum((y(1:end-1,:)+y(2:end,:)).*dx)];
end

