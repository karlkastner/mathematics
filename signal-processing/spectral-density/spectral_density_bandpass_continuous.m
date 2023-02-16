% Tue 23 Nov 17:31:17 CET 2021
% Karl Kastner, Berlin
%
%%
%% function [S_bp,Sc] = spectral_density_bandpass_continous(fx,fc,order,normalize,pp)
%%
%% output :
%% S_bp : spectral density of the bandpass filter in continuos space
%%     limit case of the discrete bandpass for dx -> 0
%% Sc   : scale factor to normalize area to 1, if noramlize = true
%%
%% input :
%% f     : frequency (abszissa)
%% fc    : central freqeuncy, location of maximum on abszissa
%% order : number of times filter is applied iteratively, not necessarily integer
%% normalize : normalize area under curve int_0^inf S(f) df = 1, if not maximum S(fc) = 1
%% pp    : powers for recombination of the lowpass filter 
function [S_bp,Sc] = spectral_density_bandpass_continuous(fx,fc,order,normalize,pp)
	if (nargin()< 3 || isempty(order))
		order = 1;
	end
	if (nargin()<4 || isempty(normalize))
		normalize = 1;
	end
	if (nargin()<5)
		pp = [];
	end
	if (isempty(pp))
		% this is identical to pp = [2,1/2,1]
		% S     = (1-S_lp).*S_lp;
		S_bp = (2*fc.*fx./(fx.^2 + fc.^2)).^2;
	else
		% lowpass density
		S_lp1 = spectral_density_lowpass_continuous(fx,fc,1,false);
		% transfer function
		T_lp1 = sqrt(S_lp1);
		% bandpass density
		S_bp     = ((1-T_lp1.^pp(1)).^pp(2).*T_lp1.^pp(3)).^2;
	end

	if (length(order)>1 || order ~= 1)
		S_bp = S_bp.^order;
	end
	
	% normalize
	switch (normalize)
		case {0}
			Sc  = 1; % do not normalize
		otherwise % analytic, for limit case df = 1/L = 0
			Sc = spectral_density_bandpass_continuous_scale(fc, order, pp);
			S_bp = Sc.*S_bp;
	end % switch
end % spectral_density_bandpass_continuous

