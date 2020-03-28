% Sat Nov 22 10:54:44 CET 2014
% Karl Kastner, Berlin

function [s obj] = std(obj)
	s = obj.stdS(obj.h,obj.edge,obj.method);
end

