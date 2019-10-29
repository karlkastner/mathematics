% Wed Jul 11 18:42:16 MSK 2012
% Karl Kästner, Berlin

function [P T BC] = promote_1d_2_3(P, T, BC)
	% BC does not change
	BC = BC;
	lp = size(P,1);
	lt = size(T,1);
	P = [P; 0.5*(P(T(:,1)) + P(T(:,2)))];
	T = [T lp+(1:lt)'];
end

