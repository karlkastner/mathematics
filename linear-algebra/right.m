% 2016-09-10 12:51:23.225201635 +0200
% Karl Kastner, Berlin
%
%% get right column by shifting columns to left
%% extrapolate or linearly exapolate rightmost column
function x=right(x,iscircular)
	if (size(x,2)<2)
		warning('right on colum vector');
	end
	if (nargin()<2 || ~iscircular)
		x = horzcat(x(:,2:end),2*x(:,end)-x(:,end-1));
	else
		x = [x(:,2:end),x(:,1)];
	end
end
