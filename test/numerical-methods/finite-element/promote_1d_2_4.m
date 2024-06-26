% Wed Jul 11 19:11:42 MSK 2012
% Karl Kästner, Berlin

function [P T BC] = promote_1d_2_4(P, T, BC)
	% BC does not change
	BC = BC;
	lt = size(T,1);
	lp = size(P,1);
	P = [P; 1/3*(2*P(T(:,1)) + P(T(:,2))); 1/3*(P(T(:,1)) + 2*P(T(:,2)))];
	T = [T lp+(1:lt)' lp+(lt+1:2*lt)'];
end

