% Wed 18 May 11:01:36 CEST 2016
% Karl Kastner, Berlin
%
%% quadratic interpolation returning value and derivative(s)
function [y0, dy_dx0] = interp1_slope(x,ys,x0)
	for idx=1:length(x0)
		[void fdx] = min((x-x0(idx)).^2);
		n = length(ys);
		l = max(1,min(fdx-1,n-2));
		r = l+2;
		dy_dx = (ys(r)-ys(l))./(x(r)-x(l));
		A = vander_1d(x(l:r)-x0(idx),2);
		c = A \ ys(l:r);
		y0(idx)       = c(1);
		dy_dx0(idx)   = c(2);
		d2y_dx20(idx) = 2*c(3);
	end
end

