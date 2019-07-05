% Di 1. Dez 13:54:19 CET 2015
% Karl Kastner, Berlin
%
%% dot product
% note : this originally gave the enclosed angle, but this is not the dot product
function dot = dot(a,b) %dx1,dy1,dx2,dy2)
	dot = sum(a.*b);
%	if (~issym(dx1))
%		cosa = (dx1.*dx2+dy1.*dy2)./(hypot(dx1,dy1).*hypot(dx2,dy2));
%	else
%		cosa = (dx1.*dx2+dy1.*dy2)./sqrt((dx1.^2+dy1.^2).*(dx2.^2+dy2.^2));
%	end
end

