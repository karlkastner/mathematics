% 2016-09-10 12:51:23.225201635 +0200
%% get right column by shifting columns to left
%% extrapolate rightmost column
function x=right(x)
	if (size(x,2)<2)
		warning('right on colum vector');
	end
	x = horzcat(x(:,2:end),2*x(:,end)-x(:,end-1));
end
