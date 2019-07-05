% Mon 22 Jan 11:08:04 CET 2018
%
%% fits coefficients of a continuous fourier transform,
%% but stores them as resampled values
%
function [vals, serr, cnt] = fourier_resampled_fit(T,t0,t,val)
	[c, serr, cnt] = fourier_fit(T, t0, t, val);
	% number of coefficients
	nc = 1+2*length(T);
	% support points for resampling
	% note : zero time is shifted to t0
	%ts   = t0(1) + (t0(end)-t0(1))*((0:no-1)/no + 0.5)';
	ts   = (t0(end)-t0(1))*((0:nc-1)/nc + 0.5)';
	F    = fourier_matrix(T,ts);
	% resample
	vals = F*c;
end

