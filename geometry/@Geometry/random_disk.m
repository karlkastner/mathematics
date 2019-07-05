% Sun 14 Jan 17:39:53 CET 2018
%
%% draw random points on the unit disk
function [P0, r, theta] = random_disk(n)
	% sqrt for equal probability, see wolfram
	r = sqrt(rand(n,1));
	theta = 2*pi*rand(n,1);

	x = r.*cos(theta);
	y = r.*sin(theta);
	P0 = [x,y];
end
