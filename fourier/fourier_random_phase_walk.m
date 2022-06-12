% 2021-08-28 14:43:20.338613448 +0200
% Karl Kästner, Berlin
%
%% evaluete fourier series where the phase undergoes a brownian motion
%
% function [x,y,phi] = fourier_random_phase_walk(a,b,f,s,L,dx)
function [x,y,phi] = fourier_random_phase_walk(a,b,f,s,L,dx,m,isbridge)
	a = cvec(a);
	b = cvec(b);
	f = rvec(f);
	if (nargin()<7)
		m = 1;
	end
	if (nargin()<8)
		isbridge = false;
	end

	% create displaced series
	nx= round(L/dx);
	x = dx*(0:nx-1)';
	%x = (0:dx:L-dx)';
	%n = length(x);
	%phi = 2*pi*1/sqrt(dx)*s(idx)*cumsum(dx*randn(length(xi),1));
	if (~isbridge)
		ex = sqrt(dx)*s*cumsum(randn(nx,m)) + 1/f*rand(1,m);
	else
		ex = sqrt(dx)*s*cumsum(randn(nx+1,m)) + 1/f*rand(1,m);
		%x_ = dx*(0:nx-1)';
		k = (0:nx)'/nx;
		% brownian bridge
		ex = ex - k.*(ex(end,:)-ex(1,:));
		ex = ex(1:nx,:);
	end

	phi = 2*pi*ex;
	if (0) %nargin()>6)
		el = randn(nx,1);
		invert = true;
		order = 2;
		el = lowpass1d_implicit(el,lpf,order,invert); 
		%el = s*el/rms(el);
		f0 = 1;
	        el = 0.275/f0*el/rms(el);
		%phi = phi + 2*pi*el;
		el = 1+el;
	end

	% ab = hypot(a,b);
	y = (   cos((x + ex).*(2*pi*f))*a ...
              + sin((x + ex).*(2*pi*f))*b ...
	    );
	% TODO compute exact minimum
	y = y-min(y);
	if (0)
		%nargin()>6)
		% must come after subtraction of mimimum
		% this shows, that this is not (only) amplitude variation,
		% modulation of the mean component,
		% amplitude of other frequencies can vary, but doe not
		% introduce the lobe at zero
		y = y.*el;
	end
end % fourier_random_phase_walk

