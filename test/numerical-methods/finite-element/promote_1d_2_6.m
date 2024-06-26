% Wed Jul 11 19:19:27 MSK 2012
% Karl Kästner, Berlin

function [P T BC] = promote_1d_2_6(P, T, BC)
	% BC does not change
	BC = BC;
	lt = size(T,1);
	lp = size(P,1);
	P = [P; 0.2*(4*P(T(:,1)) + P(T(:,2))); 0.2*(3*P(T(:,1)) + 2*P(T(:,2))); 0.2*(2*P(T(:,1)) + 3*P(T(:,2))); 0.2*(P(T(:,1)) + 4*P(T(:,2)))];
	T = [T lp+(1:lt)' lp+(lt+1:2*lt)' lp+(2*lt+1:3*lt)' lp+(3*lt+1:4*lt)'];
end % promote_1d_2_6

