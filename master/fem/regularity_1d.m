% Tue Jun 19 15:53:41 MSK 2012
% Karl Kästner, Berlin

function [h_side C] = regularity_1d(P,T,Bc)

	% lenght of element sides
	h_side = abs(P(T(:,1)) - P(T(:,2)));

	% centre coordinates
	C = 0.5*(P(T(:,1)) + P(T(:,2)));
end % regularity_2d
