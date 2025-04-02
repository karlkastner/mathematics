% Sat 26 Jun 23:13:56 CEST 2021
% Karl Kastner, Berlin
%
%% function [S_bp, S_bp1] = spectral_density_bandpass_discrete(fx,arg2,order,dx,varargin)
%% spectral density of the discrete spatial (two-sided) bandpass filter
%% in discrete space
%%
function [S_bp, S_bp1] = bandpass1d_discrete_pdf(fx,arg2,order,dx,normalize,varargin)
	if (nargin()<3 || isempty(order))
		order = 1;
	end
	if (nargin()<5 || isempty(normalize))
		normalize = true;
	end
if (0)
	[S_lp, S_lp1] = spectral_density_lp(fx,arg2,dx,1,varargin{:});
	% fourier transform of the autocorrelation function of a first order bandpass
	% scale, so that maximum is 1
	S_bp1 = 4*(1-S_lp1).*S_lp1;
else
	r      = bandpass_arg(arg2,dx,varargin{:})
	% fourier transform of the autocorrelation function of a first order bandpass
	%S_lp1  = (1 - 2*rho + rho.^2)./(1-2*rho*cos(2*pi*f)+rho.^2);
	omega  = 2*pi*fx;
	S_bp1  = 8*r.*(cos(omega*dx)-1)./(1 + 2*r.*(1-cos(omega*dx))).^2;
end

	% for the spectral density, the fourier transform has to be squred
	S_bp2 = S_bp1.^2;
	% note : the double step is important to avoid powers of negative values!
	S_bp = S_bp2.^(order);

	switch (normalize)
	case {-1}
		% nothing to do
	otherwise
		I = spectral_density_area(fx,S_bp);
		S_bp = S_bp./I;
	end
end % bandpass2d_discrete_pdf

