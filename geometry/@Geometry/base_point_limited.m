% Sun 15 Jul 13:04:33 CEST 2018
% closest point on a line section
%
%% base point (Fusspunkt) of a point on a line
function [xb,s,ds] = base_point_limited(x0,x1,x2)
	[xb, s] = Geometry.base_point(x0,x1,x2);
	ldx        = (s<0);
	xb(1,ldx)  =  x1(1);
	xb(2,ldx)  =  x1(2);
	rdx        = (s>1);
	xb(1,rdx)  =  x2(1);
	xb(2,rdx)  =  x2(2);
	
	if (nargout() > 2)
		dxy = x0 - xb;
		ds = hypot(dxy(1,:),dxy(2,:));
		%x0(:,1)-xb(1),x0(:,2)-xb(:,2));
	end
end

