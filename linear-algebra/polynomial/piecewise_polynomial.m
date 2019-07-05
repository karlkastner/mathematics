% 2015-06-09 10:03:43.428452134 +0200
%% evaluate piecewise polynomial
function y = piecewise_polynomial(x,x0,y0)
	y = zeros(size(x));
	dy0 = diff(y0);
	dx0 = diff(x0);
	dy0_dx0 = dy0./dx0;
	for idx=1:length(x0)-1
		y = y + (x >= x0(idx) & x<x0(idx+1)).*((x - x0(idx))*dy0_dx0(idx) + y0(idx));
	end
end

