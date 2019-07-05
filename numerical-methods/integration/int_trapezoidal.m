% 2018-08-29 14:27:21.089395412 +0200
%% integrate y along x with the trapezoidal rule
function int = int_trapezoidal(x,y)
	dx = diff(x);
	int = 0.5*sum((y(1:end-1)+y(2:end)).*dx);
end

