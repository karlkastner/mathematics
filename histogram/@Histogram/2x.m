% 2014-11-22 21:06:10.542874324 +0100
% Karl Kastner, Berlin
function [X obj] = 2x(obj,n)
	h = obj.h;
	N = h*n;
	X = [];
	for idx=1:length(N)
		X = [X; repmat(idx,round(N(idx)),1)];
	end
end
