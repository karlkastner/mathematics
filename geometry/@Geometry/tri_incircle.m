% Sat 10 Sep 10:53:52 CEST 2016
% Karl Kastner, Berlin
%% incircle of a triangle
function [x0, y0, R] = tri_incircle(x,y)
	% get edge lengths
	xl = left(x);
	yl = left(y);
	xr = right(x);
	yr = right(y);
	if (~issym(x))
		l  = hypot(xr-xl,yr-yl);
	else
		l = sqrt((xr-xl).^2 + (yr-yl).^2);
	end
	s  = sum(l,2);
	% 
	x0 = sum(l.*x,2)./s; % (rvec(l)*cvec(x));
	y0 = sum(l.*y,2)./s; % (rvec(l)*cvec(y));
	if (nargout() > 2)
		R = 0.5*sqrt( (l(:,2)+l(:,3)-l(:,1)) ...
			    .*(l(:,3)+l(:,1)-l(:,2)) ...
			    .*(l(:,1)+l(:,2)-l(:,3))./s );
	end
end

