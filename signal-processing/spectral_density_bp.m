%Sat 26 Jun 23:13:56 CEST 2021
% spectral density of a bandpass
% function [S_bp, S_bp1] = spectral_density_bp(f,arg2,L,n,order,varargin)
function [S_bp, S_bp1] = spectral_density_bp(f,arg2,L,n,order,varargin)
	if (nargin()<4 || isempty(order))
		order = 1;
	end
if (0)
	[S_lp, S_lp1] = spectral_density_lp(f,arg2,L,n,1,varargin{:});
	% fourier transform of the autocorrelation function of a first order bandpass
	% scale, so that maximum is 1
	S_bp1 = 4*(1-S_lp1).*S_lp1;
else
	r = filter_arg(arg2,L,n,varargin{:});

	% fourier transform of the autocorrelation function of a first order bandpass
	%S_lp1  = (1 - 2*rho + rho.^2)./(1-2*rho*cos(2*pi*f)+rho.^2);
	omega  = 2*pi*f;
	S_bp1  = 8*r*(cos(omega*L/n)-1)./(1 + 2*r*(1-cos(omega*L/n))).^2;
end

	% for the spectral density, the fourier transform has to be squred
	S_bp = S_bp1.^(2*order);
end 

