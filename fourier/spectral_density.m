% 2016-07-13 15:46:31.221013136 +0200
% Karl Kastner, Berlin
%
%% spectral density
%
% TODO explain
function [P, Y] = spectral_density(t,x,f0)
	t = cvec(t);
	x = cvec(x);
	p = zeros(size(f0));
	fdx = isfinite(x);
	t_ = t(fdx);
	x_ = x(fdx);
	% TODO 0 (mean)
	for idx=1:length(f0)
		s = sin(2*pi*t_*f0(idx));
		c = cos(2*pi*t_*f0(idx));
		% TODO make orthogonal if there are nans
		Y(idx) = 1./(c'*c)*(c'*x_) + 1i./(s'*s)*(s'*x_);
	end
	P = abs(Y);
end

