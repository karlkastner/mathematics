% Mon 22 Jan 11:08:21 CET 2018
% Karl Kastner, Berlin
%
%% interpolates a continuous fourier series that has been stored as values
%% at their support points
%
%
% TODO, this is nothing more than fourier interpolate
function [val] = fourier_resampled_predict(T,t0,vals,t)
	% number of coefficients
	nc = 1+2*length(T);
	% support points for resampling
	% note : zero time is shifted to t0 here
	%ts   = t0(1) + (t0(end)-t0(1))*((0:nc-1)/nc + 0.5)';
	ts   = (t0(end)-t0(1))*((0:nc-1)/nc + 0.5)';
	F    = fourier_matrix(T,ts);
	% get coefficients from resampled values
	c = F \ vals;

	[val] = fourier_predict(T, t0, c, t);
end % fourier_resampled_predict

