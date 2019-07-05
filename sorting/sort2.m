% 2016-01-06 13:32:12.923502801 +0100
% Karl Kastner, Berlin
%
%% sort two numbers
function [as, bs] = sort2(a,b)
	if (a <= b)
		as = a;
		bs = b;
	else
		as = b;
		bs = a;
	end
end

