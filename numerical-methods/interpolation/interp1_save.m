% 2018-02-24 13:13:27.270786778 +0100
%% make interpolation save to round off errors
%% the matlab internal interpolation suffers from rounding errors, which
%% are unacceptable when values of X and Y are large (for example UTm coordinates)
%% this normalization prevents this
function Yi = interp1_save(X,Y,Xi,varargin)
	x0  = X(1);
	X   = X-x0;
	Xi  = Xi-x0;
	L = X(end);
	X = X/L;
	Xi=Xi/L;
	
	Yi = interp1(X,Y,Xi,varargin{:});
	
end
