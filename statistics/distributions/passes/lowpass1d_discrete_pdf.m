% Sat 26 Jun 21:04:19 CEST 2021
% c.f. berkeley lecture notes, this is just the squared
%
% spectral density (S) of the spatial (two-sided) pth-order recursive lowpass
% identical power spectrum of the ar1 process of pth-recursive order
% S = |F|^p, where F is fourier transform of the autocorrelation function of a first order lowpass
%
% function [S_lp,S_lp1] = spectral_density_lowpass(f,arg2,order,dx,varargin)
function [S_lp, S_lp1] = lowpass1d_continuous_pdf(fx,arg2,order,dx,normalize,varargin)
	if (nargin()<3||isempty(order))
		order = 1;
	end
	if (nargin()<5)
		normalize = true;
	end
	r      = bandpass_arg(arg2,dx,varargin{:});
	% note that lp_arg (f,p) = bp_arg(f), i.e. the higher order is cancelled
	% by the requirement S_lp(f) = 0.5
	%r      = lowpass_arg(arg2,dx,1,varargin{:});
	omega  = 2*pi*fx;
	S_lp1  = 1./(1 + 2*r*(1-cos(omega*dx)));

	% for the spectral density, the fourier transform has to be squared
	% squaring for two-sided application
	S_lp1 = S_lp1.^2;
	% raise to power order
	S_lp  = S_lp1.^(order);
	
	switch (normalize)
	case {-1}
		% do not normalize
	otherwise
		% normalize numerically
		%df = 1/L;
		%S_bp = 2*S_bp/sum(S_bp*df);
		I = spectral_density_area(fx,S_lp);
		S_lp = S_lp./I;
	end
end 

