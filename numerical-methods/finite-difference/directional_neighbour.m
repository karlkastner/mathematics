% Sun 11 Jul 17:45:23 CEST 2021
% supports set-up of the 2D derivative matrix in a 3x3 kernel
% into an arbitrary direction
%
% this is different from p*D1x + (1-p)*D1y, which has zeros in the corners
function [r,l] = directional_neighbour(a)

% TODO it is better to just rotate the 3x3 kernel, and then to use it in the sparse matrix setup
% TODO, compare to imrotate, this is not correct
x   = ones(3,1)*[1 2 3];
y   = x';
R   = rot2(a);
xy  = R'*[flat(x),flat(y)];
xy_ = floor(xy_);
p   = xy - xy_;



	c = sin(a);
	s = sqrt(1-c*c);

	if (c>0)
	if (s>0)
	% this is only the + half
	r = [0,     (1-c)*s,   c*s;
	     0     (1-c)*(1-s),   c*(1-s);
	     0, 0, 0];
	l = [0, 0, 0;
	     c*(1-s), (1-c)*(1-s), 0;
	     c*s,     (1-c)*s, 0];
	else
		s_ = -s;
	r = [   0,               0,        0;
	        0,    (1-c)*(1-s_), c*(1-s_);
	        0,        (1-c)*s_,    c*s_];
	r = [       c*s_,        (1-c)*s_,    0;
	        c*(1-s_),    (1-c)*(1-s_),    0;
	        0,    0,    0];
	end
	else
		error('here');
	end

if (0)
	di = c+1;
	dj = s+1;
%	l = zeros(3,);
	di = ceil(di);
	dj = ceil(dj);
	l(di+2,dj+2);
end
end

