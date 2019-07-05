% Thu Nov 20 14:22:40 CET 2014
% 2015-10-30 19:18:30.973286724 +0100
% Karl Kastner, Berlin
%
%% curvature of a function in two dimensions
% 
% function [c, R] = curvature_(x,y)
% or:
% function [c, R] = curvature_(dx_ds,dy_ds,d2x_ds2,d2y_ds2)
%
%
function [c, R] = curvature(varargin)
	switch (length(varargin))
	case {2,4}
		[c,R] = curvature_(varargin{:});
	otherwise
		error('curvarture: incorrect number of input arguments');
	end	
end

% note that s is not needed, because it cancels out of the equation
function [c, R] = curvature_(dx_ds,dy_ds,d2x_ds2,d2y_ds2)
	fftflag = false;
	if (2 == nargin())
		
		if (fftflag)
		X    = cvec(dx_ds);
		Y    = cvec(dy_ds);
		n = length(X);
		D1    = spectral_derivative_1d(n,1,0);
		D2    = spectral_derivative_1d(n,2,0);
		dx    = real(ifft(D1.*fft(X)));
		dy    = real(ifft(D1.*fft(Y)));
		d2x   = real(ifft(D2.*fft(X)));
		d2y   = real(ifft(D2.*fft(Y)));

%		% resample
%		n = length(X);
%		ni = ceil(log(n/log(2)));
%		% TODO actually X(end)-h
%		Xi = linspace(X(1),X(end),ni)';
%		XYi = interp1(1:n,
		% TODO higher accuracy as a piecewise hermite polynomial?

		else
		% arg1 and arg2 are x and y not derivatives
		X    = cvec(dx_ds);
		Y    = cvec(dy_ds);
%		X = gaussfilt1(X,10);
%		Y = gaussfilt1(Y,10);

		dx   = cdiff(X);
		dy   = cdiff(Y);
		d2x  = cdiff(X,2);
		d2y  = cdiff(Y,2);
		% S drops out
		%S    = [0; cumsum(hypot(dx,dy))];
%		S    = [cumsum(hypot(dx,dy))];
%		idS  = 1./cdiff(S);
%		dx_ds   = dx.*idS;
%		dy_ds   = dy.*idS;
%		d2x_ds2 = d2x.*idS.*idS;
%		d2y_ds2 = d2y.*idS.*idS;
		end % if fftflag
	end % if 2==nargin
	% curvature
%	c =   (dx_ds.*d2y_ds2 - dy_ds.*d2x_ds2) ...
%            ./(dx_ds.*dx_ds + dy_ds.*dy_ds).^(3/2);
	c = (dx.*d2y - dy.*d2x)./(dx.*dx + dy.*dy).^(3/2);
	% radius of curvature
	R = 1./c;
end

