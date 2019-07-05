% 2011-04-14 20:55:52.000000000 +0200 geometry/ellipseX.m
%
%% x-coordinates of y-coordinates of an ellipse
function [X, Y] = ellipseX(a,b,c,d,e,f,Y1)
	X1 = -(b*Y1 - d + (b^2*Y1.^2 - 2*b*d*Y1 + d^2 - 4*a*c*Y1.^2 + 4*a*e*Y1 - 4*a*f).^(1/2))/(2*a);
	X2 =  (d - b*Y1 + (b^2*Y1.^2 - 2*b*d*Y1 + d^2 - 4*a*c*Y1.^2 + 4*a*e*Y1 - 4*a*f).^(1/2))/(2*a);
	if (size(X1,1) > size(X1,2))
		Y = [Y1; Y1];
		X = [X1; X2];
	else
		Y = [Y1 Y1];
		X = [X1 X2];
	end
end % ellipseY

