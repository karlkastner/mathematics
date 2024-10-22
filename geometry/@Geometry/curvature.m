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
		error('Curvarture:IncorrectNumberOfInputArguments","Incorrect number of input arguments');
	end	
end

% note that s is not needed, because it cancels out of the equation
function [c, R] = curvature_(dx_ds,dy_ds,d2x_ds2,d2y_ds2)
	fftflag = false;
	if (2 == nargin())
		if (fftflag)
			% Note, this works only well if ds is smoothly varying
			X    = cvec(dx_ds);
			Y    = cvec(dy_ds);
			n     = length(X);
			D1    = spectral_derivative_1d(n,1,0);
			D2    = spectral_derivative_1d(n,2,0);
			dx    = real(ifft(D1.*fft(X)));
			dy    = real(ifft(D1.*fft(Y)));
			d2x   = real(ifft(D2.*fft(X)));
			d2y   = real(ifft(D2.*fft(Y)));
		else
			% arg1 and arg2 are x and y not derivatives
			X    = cvec(dx_ds);
			Y    = cvec(dy_ds);
			dx   = cdiff(X);
			dy   = cdiff(Y);
			d2x  = cdiff(X,2);
			d2y  = cdiff(Y,2);
		end % if fftflag
	end % if 2==nargin
	% curvature
	c = (dx.*d2y - dy.*d2x)./(dx.*dx + dy.*dy).^(3/2);
	% radius of curvature
	R = 1./c;
end

