% Di 2. Feb 18:22:17 CET 2016
% Karl Kastner, Berlin
%% detrending by polynomial regression
function [Y slope] = detrend(X,Y,W)
	if (nargin() < 3)
		W = [];
	end
	slope = PolyOLS.slope(X,Y,W);
	Y = Y - slope*X;
end % detrend

