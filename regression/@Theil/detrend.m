% Di 2. Feb 18:22:17 CET 2016
% Karl Kastner, Berlin
%% linear detrending of a set of samples by the Theil-Senn Slope
function [Y, slope] = detrend(obj,X,Y,W)
	if (nargin() < 3)
		W = [];
	end
	slope = obj.slope(X,Y,W);
	Y = Y - slope*X;
end % detrend

