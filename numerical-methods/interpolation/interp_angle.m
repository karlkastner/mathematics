% Mon 23 Jan 17:19:55 CET 2017
% Karl Kastner, Berlin
%
%% interpolate an angle
%function ai = interp_angle(t,a,ti,varargin)
function ai = interp_angle(t,a,ti,varargin)
	s  = sin(a);
	c  = cos(a);
	si = interp1(t,s,ti,varargin{:});
	ci = interp1(t,c,ti,varargin{:});
	ai = atan2(ci,si);
end

