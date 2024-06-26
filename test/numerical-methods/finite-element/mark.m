% Tue Jun 26 19:57:37 MSK 2012
% Karl Kästner, Berlin

function M = mark(err, thresh)
	% mark elements for refinement
	M = find(err > thresh);
end % mark_2d

