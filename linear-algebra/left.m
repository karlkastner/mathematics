% 2016-09-10 12:51:21.485248654 +0200
% Karl Kastner, Berlin
%
%% left element of vector or left column of matrix
%% identical to shifting right
%% leftmost element is extrapolated or circulated
function x = left(x,iscircular)
	if (size(x,2)<2)
		warning('left on colum vector');
	end
	if (nargin()<2 || ~iscircular)
		x = horzcat(2*x(:,1)-x(:,2),x(:,1:end-1));
	else
		x = [x(:,end), x(:,1:end-1)];
	end
	%x = horzcat(x(:,end),x(:,1:end-1));
end

