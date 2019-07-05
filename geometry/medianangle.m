% 2016-07-09 12:49:22.485242539 +0200
%
%% median angle
%% angle, that has the smallest squared distance to all others
%
function me = median_angle(a)
	me = a(1);
	s = inf(size(a));
	for idx=1:length(a)
		d       = a-a(idx);
		d       = wrapToPi(d);
%		flag    = d < -pi;
%		d(flag) = 2*pi + d(flag);
%		flag    = d > pi;
%		d(flag) = 2*pi-d(flag);
		s(idx) = abs(sum(d>0)-sum(d<0));
	end
	% TODO, treat ties
	[s sdx] = min(s);
	me = a(sdx);
end

