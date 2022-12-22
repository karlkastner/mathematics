% Mon 19 Dec 17:03:02 CET 2022
% gaussian window in 2D, allowing for non-integer radii
% TODO allow for anisotropy
function [w, dof] = gausswin2(n, nf)
	[fx, fy, fr, ft] = fourier_axis_2d([1,1],n);

	% Gaussian window
	w = normpdf(fr,0,nf/log(4));

	% normalize sum of weigths to 1
	w = w/sum(w(:));

	% degrees of freedom of window
	% when values are uncorrelated
	dof = sum(w(:)).^2./(sum(w(:).^2));
end

