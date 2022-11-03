% Fri 22 Apr 13:28:53 CEST 2022
%
% function S = spectral_density_bandpass_2d(fr,fc,p,normalize,q)
%
% fr : radial frequency sqrt(fx^2 + fy^2)
% Lr : cut-off distance
function [S_bp,Sc] = spectral_density_bandpass_2d(fr,fc,p,normalize,pp)
	if (nargin()<4)
		normalize = 1;
	end
	if (nargin()<5)
		pp = [];
	end
	%fc = fc/(2*2^(1/3) - 1)^(1/2);
	if (isempty(pp))
		S_bp = (2*fc*fr./(fr.^2 + fc.^2)).^2;
		%S = (4*(1 - (1 + (fr./fc).^2).^(3/4))).^2./(1 + (fr./fc).^2).^3;
	else
		%S = (4*(1 - ((fr./fc).^2 + 1).^(3/2)).^2)./(1 + (fr/fc).^2).^6;
		% lowpass density
		S_lp1 = spectral_density_lowpass_continuous(fr,fc,1,false);
		% transfer function
		T_lp1 = sqrt(S_lp1);
		% bandpass density
		S_bp     = ((1-T_lp1.^pp(1)).^pp(2).*T_lp1.^pp(3)).^2;
	end

	if (nargin()>2 && ~isempty(p))
		S_bp = S_bp.^p;
	end
	switch(normalize)
		case {0} % nothing
			Sc = 1;
		otherwise
			Sc = spectral_density_bandpass_2d_scale(fc,p);
			S_bp = Sc.*S_bp;
	end % switch
end % spectral_density_bandpass_2d

