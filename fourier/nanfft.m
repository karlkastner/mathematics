% 2016-07-13 15:29:54.588967826 +0200
%
%% discrete fourier transform of a data series with gaps
%
function [f, Fi] = nanfft(x,m)
	n = length(x);
	% get inverse Fourier matrix
	% TODO one could get directly the qr factors here
	Fi = idftmtx_man(n,m);
	% determine valid samples
	fdx = isfinite(x);
	% determine coefficients with least squares
	f = Fi(fdx,:) \ x(fdx);
end

