% Tue 23 Jan 13:19:22 CET 2018
% Karl Kastner, Berlin
%
%% rotation matrix from direction vector
% R   = rot2dir(dir) = rot2(atan2(dir(2),dir(1)))
% dir = R*[1;0]
function R = rot2dir(dir)
	dir = dir/norm(dir);
	c = dir(1);
	s = dir(2);
	R =  [ c, -s;
               s, c];
end % rot2dir

