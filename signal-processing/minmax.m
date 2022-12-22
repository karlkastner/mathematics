% Tue 13 Dec 15:43:42 CET 2022
function [mindx,maxdx] = minmax(x)
	maxdx_=1;
	mindx_=1;
	mindx = 1;
	maxdx = 1;
	% find min
	for idx=1:length(x)
		if (x(idx) < x(mindx_))
			mindx_ = idx;
			maxdx_ = idx;
		elseif (x(idx) > x(maxdx_))
			maxdx_ = idx;
			if ((x(maxdx_)-x(mindx_)) > (x(maxdx)-x(mindx)))
				maxdx = maxdx_;
				mindx = mindx_;
			end
		end
	end
end

