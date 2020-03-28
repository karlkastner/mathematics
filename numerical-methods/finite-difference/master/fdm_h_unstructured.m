% Thu May 24 19:38:14 MSK 2012
% Karl Kästner, Berlin

function [Hxl Hxr Hyl Hyr] = fdm_h_unstructured(P)
	% left
	Hxl = abs(P(idx,1) - P(N(idx,1),1));
	% right
	Hxr = abs(P(N(idx,2),1) - P(idx,1));
	% bottom
	Hyl = abs(P(idx,2) - P(N(idx,3),2));
	% top
	Hyr = abs(P(N(idx,4),1) - P(idx,2));
end % fdm_h_unstructured()

