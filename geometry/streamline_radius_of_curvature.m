% Sat 10 Mar 14:37:35 CET 2018
%
%% streamline radius of curvature
%% simplifies when rotatate to streamwise coordinates to R = 1/dv/ds * u
% function R = streamline_radius_of_curvature(u,du_dx,du_dy,v,dv_dx,dv_dy)
%
% note that this equation fails when u==v, also for kalkwijk
function R = streamline_radius_of_curvature(u,du_dx,du_dy,v,dv_dx,dv_dy)
	% abad 2007 (typo in their denominator was corrected)
	% also see derivation script
	R =    (u.^2 + v.^2).^(3/2) ...
	    ./ (dv_dx.*u.^2 + (-du_dx + dv_dy).*u.*v - du_dy.*v.^2);
end

