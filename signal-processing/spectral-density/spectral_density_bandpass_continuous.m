% Tue 23 Nov 17:31:17 CET 2021
% Karl Kastner, Berlin
%
%% S : spectral density of the bandpass filter in continuos space
%%     limit case of the discrete bandpass for dx -> 0
%%
%% function [S,IS] = spectral_density_bandpass_continous(fx,fc,order)
%% f     : frequency (abszissa)
%% fc    : central frequncy, location of maximum on abszissa
%% order : number of times filter is applied iteratively, not necessarily integer
function [S,IS] = spectral_density_bandpass_continuous(fx,fc,order,normalize,q)
	if (nargin()< 3||isempty(order))
		order = 1;
	end
	if (nargin()<4 || isempty(normalize))
		normalize = 1;
	end
	if (nargin()<5)
		q = 1;
	end
	switch (q)
	case {1}
		% from recombining a first order lowpass
		% this is the spectral density scaled so that the max at f=fc is 1
		S = (4*(1 - sqrt(1 + 3*(fx./fc).^2))./(1 + 3*(fx./fc).^2)).^2;
	case {2}
		% from recombining a second order lowpass
		% this expression is simpler
		S = (2*fc*fx./(fx.^2 + fc.^2)).^4;
	otherwise
		error('not yet implemented');
	end

	% nb: power 4 has to be separated, to ensure positivity
	if (order ~= 1)
		S = S.^order;
	end
	
	% normalize
	switch (normalize)
	case {-1}
		% do not normalize
		if (issym(S))
			syms Sc
		else
			Sc = 1;
		end
	case {0} % analytic, for limit case df = 1/L = 0
		%IS = spectral_density_bandpass_continuous_scale(fc, order);
		%Sc = spectral_density_bandpass_continuous_scale(fc, order, 0,1);
		Sc = spectral_density_bandpass_continuous_scale(fc, order, 1,q)
		S = Sc.*S;
	otherwise % numerical, for finite L
		%Sc = 1./spectral_density_area(fx, S);
		Sc = spectral_density_bandpass_continuous_scale(fc, order, 1, q)
%		L = 100;
%		Sc = 1./quad(@(fx) spectral_density_bandpass_continuous(fx,fc,order,-1),0,L*fc);
%		Sc = spectral_density_bandpass_continuous_scale(fc, order);
%		df = fc/100;
%		1./sum(df*spectral_density_bandpass_continuous(0:df:(L*fc),fc,order,-1))

	end
	S = Sc.*S;
end % spectral_density_bandpass_continuous

