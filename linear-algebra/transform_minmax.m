% Tue 24 Dec 09:44:58 +08 2019
function x = transform_min_max(x,xmin_,xmax_)
	xmin = min(x);
	xmax = max(x);
	if (xmin == xmax)
		x = 0.5*(xmin_ + xmax_);
	else
		x = xmin_ + ((xmax_-xmin_)/(xmax-xmin))*(x-xmin);
	end
end

