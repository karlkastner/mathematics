% Thu Jun 14 16:45:57 MSK 2012
% Karl Kästner, Berlin

function [area h_side s_angle C] = recalculate_regularity_2d(selection, P, T, area, h_side, s_angle, C)
	for tdx=selection
		% calculate element centre coordinates
		C(tdx,:) = 1/3*(P(T(tdx,1),:) + P(T(tdx,2),:) + P(T(tdx,3),:));

		% calculate element area
                A = [   1 P(T(tdx,1),:);
                        1 P(T(tdx,2),:);
                        1 P(T(tdx,3),:) ];
		area(tdx) = 0.5*abs(det(A));

		% todo: warn for clockwise (area<0) and degenerated triangles (area/h_max < eps)
	
		% length of side opposit the points
		% todo, vectorise
		a = sqrt((P(T(tdx,2),1) - P(T(tdx,3),1))^2 + (P(T(tdx,2),2) - P(T(tdx,3),2))^2);
		b = sqrt((P(T(tdx,1),1) - P(T(tdx,3),1))^2 + (P(T(tdx,1),2) - P(T(tdx,3),2))^2);
		c = sqrt((P(T(tdx,1),1) - P(T(tdx,2),1))^2 + (P(T(tdx,1),2) - P(T(tdx,2),2))^2);
		h_side(tdx,:) = [a b c];

		% sine of interior angles at the points
		% todo, vectorise
		sin_a = area(tdx)/(0.5*b*c);
		sin_b = area(tdx)/(0.5*a*c);
		sin_c = area(tdx)/(0.5*a*b);
		s_angle(tdx,:) = [sin_a sin_b sin_c];
	end
end % recalculate_regularity_2d

