% Mon  5 Oct 13:40:26 +08 2020
% TODO this assumes there are n^2-1 elements
% note : when np==1 (n == 3) the end points cannot be extrapolated
function x = mg_interpolate(xp)
	n            = 2*length(xp)+1;
	if (n == 3)
		x = xp*ones(3,1);
	else
	x            = zeros(n,1);
	x(1)         = 1.5*xp(1)-0.5*xp(2);
	x(2:2:end-1) = xp;
	x(3:2:end-2) = 0.5*(xp(1:end-1)+xp(2:end));
	x(end)       = 1.5*xp(end)-0.5*xp(end-1);
	end
end

