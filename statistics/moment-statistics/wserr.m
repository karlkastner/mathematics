% 2016-03-07 16:40:22.339348733 +0100
%% weighted root mean square error
function [serr] = wserr(w,x,np)
	if (nargin < 3 || isempty(np))
		np = 1;
	end
	n = wdof(w);
	serr = sqrt(wvar(w,x)/(n-np));
end % wserr

