
function  [s d c a] = check_area_2d(P,T,B)
	s = 0;
	d = 0;
	c = 0;
	a = 0;
	for idx=1:size(T,1)
		A = [	1 P(T(idx,1),:)
			1 P(T(idx,2),:)
			1 P(T(idx,3),:) ];
		area = det(A);
		if (abs(area) < 1e-12)
			c = c+1;
		end
		if (area < 0)
			a = a + 0;
		end
		s = s + abs(area);
	end
	if ( c > 0)
		'warning colinear (degenerated) triangles'
	end
	if ( a > 0 )
		'warning: clockwise triangles'
	end
	% length of boundary
	for idx=1:size(B,1)
		d = d + norm(P(B(idx,1),:) - P(B(idx,2),:));
	end
	s = 0.5*s;
end % check_area_2d

