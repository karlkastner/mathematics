% Wed Jul 11 19:14:03 MSK 2012
% Karl Kästner, Berlin

function [P T BC] = promote_1d_2_5(P, T, BC)
	% BC does not change
	BC = BC;
	lt = size(T,1);
	lp = size(P,1);
	P = [P; 0.25*(3*P(T(:,1)) + P(T(:,2))); 0.5*(P(T(:,1)) + P(T(:,2))); 0.25*(P(T(:,1)) + 3*P(T(:,2)))];
	T = [T lp+(1:lt)' lp+(lt+1:2*lt)' lp+(2*lt+1:3*lt)' ];
end % promote_1d_2_5

