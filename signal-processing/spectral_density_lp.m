% Sat 26 Jun 21:04:19 CEST 2021
% c.f. berkeley lecture notes, this is just the squared
% power spectrum of the ar1 process
function [S_lp,S_lp1] = spectral_densitylp(f,arg2,L,n,order,varargin)
	if (nargin()<4||isempty(order))
		order = 1;
	end

	r = filter_arg(arg2,L,n,varargin{:});

	% fourier transform of the autocorrelation function of a first order bandpass
	%S_lp1  = (1 - 2*rho + rho.^2)./(1-2*rho*cos(2*pi*f)+rho.^2);
	omega  = 2*pi*f;
	S_lp1  = 1./(1 + 2*r*(1-cos(omega*L/n)));
	

%	if (normalize)
%		f1 = (1-2*f0+f0^2)/
%	end

%	f1 = (1-2*f0+f0^2)./(1-2*f0*cos(pi*f)+f0^2);
	%if (nargin()>2 && ~isempty(order) && order ~= 1)

	% for the spectral density, the fourier transform has to be squared
	S_lp = S_lp1.^(2*order);
	%else
%		fo = f1;
%	end
end 


