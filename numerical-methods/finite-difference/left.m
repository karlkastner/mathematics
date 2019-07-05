% 2016-09-10 12:51:21.485248654 +0200
%% left element of vector, leftmost column is extrapolated
function x=left(x)
	if (size(x,2)<2)
		warning('left on colum vector');
	end
	x = horzcat(2*x(:,1)-x(:,2),x(:,1:end-1));
	%x = horzcat(x(:,end),x(:,1:end-1));
end

