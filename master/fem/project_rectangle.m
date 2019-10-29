% Thu Jun 14 21:38:05 MSK 2012
% Karl Kätner, Berlin

% projects the point number pdx onto the box side with dimension L0
function P = project_rectangle(pdx,P,bdx,L0,x0)
	left   = L0(1)-x0(1);
	right  = L0(1)+x0(1);
	bottom = L0(2)-x0(2);
	top    = L0(2)+x0(2);
	switch (bdx)
		% sides
		case {1}
			P(pdx,1) = left;
		case {2}
			P(pdx,1) = right;
		case {3}
			P(pdx,2) = bottom;
		case {4}
			P(pdx,2) = top;
		% corners
		case {5}
			P(pdx,:) = [left bottom];
		case {6}
			P(pdx,:) = [left top];
		case {7}
			P(pdx,:) = [right bottom];
		case {8}
			P(pdx,:) = [right top];
		otherwise
			error('','unknown boundary class');
	end	
end % project_rectangle

