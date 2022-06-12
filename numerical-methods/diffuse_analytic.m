% Mon  2 May 14:18:38 CEST 2022
% Mon 30 May 10:14:45 CEST 2022
% analytic solution to the heat equation
% dy/dt = e*D^2*y
%
% fft    : (analytic) solution via fourier transform
% no fft : analytic solution by convolution
% note   : for non-smooth solutions, variant the no-fft assures positivity
%          the "no-fft" variant is also implimented via fourier transform
function zt = diffuse_analytic(t,z0,L,e,usefft)
	% fundamental solution:
	% 1/(4*pi*t)^(n/2)*exp(-r.^2/(4*t))

%	if (nargin()>3&&~isempty(e))
%		t = e*t;
%	end

	if (nargin()<4)
		usefft = true;
	end

	if (isvector(z0))
		ndim = 1;
	else
		ndim = 2;
	end

	% convolve with Gaussian in fourier space
	switch (ndim)
	case {1}
		if (usefft)
			or = 2*pi*fourier_axis(L,length(z0));
			or = cvec(or);
			fz0 = fft(z0);
			fg = exp(-(or.*or)*t);
			fz = fg.*fz0;
			zt = ifft(fz);
		else
			zt = diffuse(e*t,z0,L);
		end
	case {2}
		if (usefft)
			n = size(z0);
			[fx, fy, fr, ft] = fourier_axis_2d(L,n);
			or = 2*pi*fr;
			dx = rvec(L)./rvec(n);
			fz = fft2(z0);
			fg = exp(-(or.*or).*t);
			fz = fg.*fz;
			zt = ifft2(fz);
			%s = sum(zt(:));
			%zt = zt./s;
		else
			% the kernel is gaussian, so the dimensions are separable
			% convolve along x
			zt = diffuse(e(1)*t,z0,L(1));
			% convolve along y
			zt = diffuse(e(2)*t,zt.',L(2)).';
		end
	otherwise
		error('not yet implemented');
	end
end

function zt = diffuse(t,z0,L)
	nx = size(z0,1);
	dx = L/nx;
	if (0==mod(length(z0),2))
		x = dx*[0:nx/2-1,-nx/2:-1]';
	else
		x = dx*[0:(nx-3)/2,-(nx+1)/2:-1]';
	end
	if (0)
		g = dx./sqrt(4*pi*t)*exp(-(x.^2)./(4*t));
	else
		g = dx/2*(erf((x+dx/2)/sqrt(4*t)) - erf((x-dx/2)/sqrt(4*t)));
%	g = dx./sqrt(4*pi*t)*exp(-((x-0.5*dx).^2)./(4*t));
%	g = g + dx./sqrt(4*pi*t)*exp(-((x+0.5*dx).^2)./(4*t));
%	g = g/2;
	end
	mg = max(g);
	if (mg>1)
		warning('max(g)>1');
	end
	g = g/sum(g);
%
%t
%clf
%t
%plot(g)
%[z0(:,1),zt(:,1)])
%pause
	if (0)
		zt = cconv(z0,g,length(z0));
	else
		zt = ifft(fft(g).*fft(z0));
	end			
end % diffuse

