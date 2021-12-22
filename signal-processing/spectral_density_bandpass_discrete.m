%Sat 26 Jun 23:13:56 CEST 2021
% spectral density of a bandpass
% function [S_bp, S_bp1] = spectral_density_bp(f,arg2,dx,order,varargin)
function [S_bp, S_bp1] = spectral_density_bandpass_discrete(fx,arg2,order,dx,varargin)
	if (nargin()<4 || isempty(order))
		order = 1;
	end
	if (length(varargin)>1)
		normalize = varargin{2};
	else
		normalize = true;
	end
if (0)
	[S_lp, S_lp1] = spectral_density_lp(fx,arg2,dx,1,varargin{:});
	% fourier transform of the autocorrelation function of a first order bandpass
	% scale, so that maximum is 1
	S_bp1 = 4*(1-S_lp1).*S_lp1;
else
	r      = filter_arg(arg2,dx,varargin{:});
	% fourier transform of the autocorrelation function of a first order bandpass
	%S_lp1  = (1 - 2*rho + rho.^2)./(1-2*rho*cos(2*pi*f)+rho.^2);
	omega  = 2*pi*fx;
	S_bp1  = 8*r.*(cos(omega*dx)-1)./(1 + 2*r.*(1-cos(omega*dx))).^2;
end

	% for the spectral density, the fourier transform has to be squred
	S_bp = S_bp1.^(2*order);

	if (~issym(fx) && normalize)
	%length(varargin)>1 && varargin{2})
		% normalize
		%df = 1/L;
		%S_bp = 2*S_bp/sum(S_bp*df);
		I = spectral_density_area(fx,S_bp);
		S_bp = S_bp./I;
	end
end % spectral_density_bp 

