% Di 1. Dez 14:05:10 CET 2015
% Karl Kastner, Berlin
%
%% cross product in two dimensions
function sina = cross2(dx1,dy1,dx2,dy2)
	sina = (dx1*dy2 - dy1.*dx2)./(hypot(dx1,dx2).*hypot(dy1,dy2));
end

