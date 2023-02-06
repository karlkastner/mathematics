% Mon 19 Dec 17:03:02 CET 2022
% circular window in 2D, allowing for non-integer radii
% TODO allow for anisotropy
function [w, nf2] = circwin(n, nf)
	[fx, fy, fr, ft] = fourier_axis_2d([1,1],n);

	% integer part
	nfi = floor(nf);
	p   = nf-nfi;

	% circular window
	w = (1-p)*(fr <= nfi) + p*(fr <= nfi+1);

	% area of window
	nf2 = sum(w(:));

	% normalize sum of weights to 1
	w = w/nf2;
end

