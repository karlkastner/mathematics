% 2015-07-19 15:46:14.176010177 +0200
%% detect peaks in a vector
%% requires function value to fall to p*max before new value is allowed
% TODO symmetry require, accept when x(i) < p/2*(x_m1 + x_m2), m1 < i < m2
% TODO, the current implementation may miss the first extremum
function fdx = detect_extreme(x,p)
	id = 1;
	s  = 1;
	fdx=[];
	for idx=2:length(x)
		if (s>0)
			if ( x(idx) > x(id) )
				id = idx;
			elseif (x(idx) < -p*x(id))
			%if (s*x(idx) < -p*s*x(id))
				% remember extremum
				fdx(end+1) = id;
				id=idx;
				s = -s;
			end
		else
	%		[idx x(idx) x(id)]
			if ( x(idx) < x(id) )
				id = idx;
			elseif (x(idx) > -p*x(id))
				% remember extremum
				fdx(end+1) = id;
				id=idx;
				s = -s;
			end
		end
	end
	% take over last extremum
	fdx(end+1) = id;
end

