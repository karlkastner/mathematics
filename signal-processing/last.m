% Thu 20 Jul 10:41:57 CEST 2017
%% lake tail, but for matrices
function [x] = last(x)
	if (isvector(x))
		x = x(end);
	else
		x = x(end,:);
	end
end
