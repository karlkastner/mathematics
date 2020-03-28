% Sat Nov 22 10:54:25 CET 2014
% Karl Kastner, Berlin

function [s2 obj] = var(obj)
	s2 = obj.varS(obj.h,obj.edge,obj.method);
%	if (strcmp(obj.method,'trapezoidal'))
%		s2 = s2 - 1/3*(edge(2)-edge(1)).^2;
%	end
end

