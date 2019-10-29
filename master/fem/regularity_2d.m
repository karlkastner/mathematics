% Wed May  2 17:45:18 MSK 2012
% Karl Kästner, Berlin

% derived triangle and boundary data
%
function [a_sum l_sum area l_boundary h_side s_angle C] = regularity_2d(P,T,Bc)
	% get length
	lt = size(T,1);

	area    = zeros(lt,1);
	h_side  = zeros(lt,3);
	s_angle = zeros(lt,3);
	C       = zeros(lt,2);

	% calculate derived triangle properties
	[area h_side s_angle C] = recalculate_regularity_2d(1:lt,P,T,area,h_side,s_angle,C);

	% boundary length
	l_boundary = sqrt( (P(Bc(:,1),1) - P(Bc(:,2),1)).^2 + (P(Bc(:,1),2) - P(Bc(:,2),2)).^2 );

	a_sum = sum(area);
	l_sum = sum(l_boundary);
end % regularity_2d

