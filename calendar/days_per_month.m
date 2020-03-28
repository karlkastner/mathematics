% 2017-06-15 17:43:19.789343497 +0200
% TODO account for leap years
function [nd d_start] = days_per_month()
	%     J  F  M A  M   J   J  A S   O   N    D
	nd = [31 28 31 30 31 30 31 31 30 31 30 31];
	d_start = cumsum([0 nd]);

end

