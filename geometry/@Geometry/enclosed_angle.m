% Wed 18 Jul 18:46:03 CEST 2018
%
%% angle enclosed between two lines
function [cosa,alpha] = enclosed_angle(a,b)
	la = hypot(a(1,:),a(2,:));
	lb = hypot(b(1,:),b(2,:));
	cosa = Geometry.dot(a,b)./(la.*lb);
	if (nargout() > 1)
		alpha = acos(cosa);
	end
end

