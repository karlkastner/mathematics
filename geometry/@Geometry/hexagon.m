% Fri 15 Jul 16:19:20 CEST 2016
% Karl Kastner, Berlin
%% coordinates of a hexagon, scaled and rotated
% xy0   : centre coordinate
% scale : diameter = twice the distance from centre to side mid point
% side length and distance from centre to vertex: scale*sqrt(3)/3
function xy = hexagon(xy0,scale,alpha)
	% 1/(2*cos(30deg)) = 1/sqrt(3)
	% sin(3) = 0.5
	xy = [1/sqrt(3)*[1, 0.5, -0.5, -1, -0.5, 0.5 ];
                    0.5*[0,    1,   1,  0,  -1,    -1]];
	% rotate
	if (nargin() > 2)
		R = rot2(alpha);
		xy = R*xy;
	end
	% scale
	if (nargin() > 1)
		xy = scale*xy;
	end
	% translate
	if (nargin()>0)
		xy(1,:) = xy(1,:)+xy0(1);
		xy(2,:) = xy(2,:)+xy0(2);
	end
end % hexagon

