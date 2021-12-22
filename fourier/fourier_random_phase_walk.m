% 2021-08-28 14:43:20.338613448 +0200
% function [x,y,phi] = fourier_random_phase_walk(a,b,f,s,L,dx)
function [x,y,phi] = fourier_random_phase_walk(a,b,f,s,L,dx)
	a = cvec(a);
	b = cvec(b);
	f = rvec(f);
	% create displaced series
	x = (0:dx:L-dx)';
	n = length(x);
	%phi = 2*pi*1/sqrt(dx)*s(idx)*cumsum(dx*randn(length(xi),1));
	ex = sqrt(dx)*s*cumsum(randn(n,1));
	phi = 2*pi*ex;
	% ab = hypot(a,b);
	y = (   cos((x + ex)*(2*pi*f))*a ...
              + sin((x + ex)*(2*pi*f))*b ...
	    );
	% TODO compute exact minimum
	y = y-min(y);
end % fourier_random_phase_walk

