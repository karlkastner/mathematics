% 2013-10-25 17:21:07 +0700
% Karl Kastner, Berlin
%
function [mx, my] = range_mean(range,x,y)
	for idx=1:size(range,1)
		fdx = find( (x >= range(idx,1)) .* (x <= range(idx,2)));
		my(idx,:) = mean(y(fdx,:),1);
		mx(idx,1) = mean(x(fdx));
	end
end

