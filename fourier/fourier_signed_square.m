% 2017-04-07 14:31:30.015413021 +0200
% Karl Kastner, Berlin
%
%% coefficients of the fourier series of | cos a + cos t | (cos a + cos t)
%% in general
%%            cos a is midrange
%%            cos t is tidal variation
%% c.f Dronkers
% TODO this does not seem to work for alpha < 0
function an = an(alpha,n)
	if (issym(alpha))
		syms pi_
	else
		pi_ = pi;
	end
	switch (n)
	case {0}
		an = 1/pi_*( (2+cos(2*alpha))*(1/2*pi_-alpha) + 3/2*sin(2*alpha));
	%	an = an/2;
	case {1}
		an = 1/pi_*(4*cos(alpha)*(1/2*pi_-alpha) + 3*sin(alpha) + 1/3*sin(3*alpha));
		%an = 1/pi_*(4*cos(alpha*(1/2*pi_-alpha)) + 3*sin(alpha) + 1/3*sin(3*alpha));
	case {2}
		an = 1/pi_*(1/2*pi_ - alpha + 2/3*sin(2*alpha) - 1/12*sin(4*alpha));
	otherwise % note: this does not apply to 1..2, division by zero!
	an = (-1)^(n+1)*2/(n*pi_)...
		*( sin((n-2)*alpha)/((n-1)*(n-2)) ...
		   - 2*sin(n*alpha)/((n+1)*(n-1)) ...
		   + sin((n+2)*alpha)/((n+2)*(n+1)) );
	end
end

