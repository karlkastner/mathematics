% 2018-10-06 16:43:42.872188172 +0200
%% segment end point to segment mid point transformation for regular 1d grids
function x_ = pwmid(x)
	n = floor(size(x,1)/2);
	x_ = 0.5*(x(1:2:2*n-1,:)+x(2:2:2*n,:));
	if (2*n<length(x))
		x_(end+1,:) = x(end,:);
	end
end

