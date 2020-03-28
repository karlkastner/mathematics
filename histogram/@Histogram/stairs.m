% Fr 15. Jan 19:39:45 CET 2016
% prepair histogram and edges for staircase plotting
% accounts for changing stap width
function [wh e obj] = stairs(obj)
	[wh e] = obj.stairsS(obj.h,obj.edge);
end

