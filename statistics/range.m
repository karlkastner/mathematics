% 2015-06-09 21:17:23.351522837 +0200
% Karl Kastner, Berlin
%
%% mid range
function [range, midrange] = range(X)
	maX = max(X);
	miX = min(X);
	range = maX - miX;
	midrange = 1/2*(maX + miX);
end

