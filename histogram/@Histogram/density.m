% 2016-03-15 20:02:49.305194853 +0100
% density for plotting at bin centres
function [density_wh, centre] = density(h,edge_A)
	n = size(h,1);
	m = size(h,2);
	% width normalisation
	w = abs(rvec(diff(edge_A)));

	centre = mid(edge_A);

	if (isvector(h))
	h=rvec(h);
	end
	w = rvec(w);

	% density
	density_wh = (h./w);

end

