% 2011-04-14 20:31:17.000000000 +0200 geometry/ellipseY.m
% y-coordinates of x-coordinates of an ellipse
function [X Y] = ellipseY(a,b,c,d,e,f,X)
	Y1 = -(b*X - e + sqrt(b^2*X.^2 - 2*b*e*X + e^2 - 4*a*c*X.^2 + 4*c*d*X - 4*c*f))/(2*c);
	Y2 =  (e - b*X + sqrt(b^2*X.^2 - 2*b*e*X + e^2 - 4*a*c*X.^2 + 4*c*d*X - 4*c*f))/(2*c);
	if (size(X,1) > size(X,2))
		Y = [Y1; Y2];
		X = [X; X];
	else
		Y = [Y1 Y2];
		X = [X X];
	end
end % ellipseY

