% Fri 14 Oct 17:11:27 CEST 2016
%
%% arc length of a two dimensional function
% TODO use arc-length
function s = integrate_path2(x,y)
	n = length(x);
	t1 = (0:(n-1))/(n-1);
	
	s1 = cumsum(hypot(diff(x),diff(y)));
if (1)
	t2 = (0:(2*n-2))/(2*n-2);
	x2 = interp1(t1,x,t2,'pspline');
	y2 = interp1(t1,y,t2,'pspline');
	s2 = cumsum(hypot(diff(x2),diff(y2)));
	s = 2*s2(2:2:end) - s1;
else
	x2 = x(1:2:end);
	y2 = y(1:2:end);
	s2 = cumsum(hypot(diff(x2),diff(y2)));
	s = 2*s1(2:2:end) - s2;	
end
%pause
end

