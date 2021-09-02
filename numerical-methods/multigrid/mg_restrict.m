% Mon  5 Oct 13:40:42 +08 2020
% TODO this assumes there are n^2-1 elements
function xp = mg_restrict(x)
	xp = 0.25*(x(1:2:end-1) + 2*x(2:2:end) + x(3:2:end));
end

