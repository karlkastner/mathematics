% Fri  7 Jan 16:11:57 CET 2022
%
%% transform mode (maxima) of the bandpass spectral density into the paramter
%% of the underlying distribution 
%
% function [p] = spectral_density_bandpass_max2par(fc,Sc,p0)
function [p] = spectral_density_bandpass_continuous_max2par(fc,Sc,p0,pp)
	if (nargin()<3)
		p0 = 1;
	end
	if (nargin()<4)
		pp = [];
	end
%	p  = fzero(@(p) spectral_density_bandpass_continuous_scale(fc,abs(p),pp) - Sc, p0);
	% n.b: lsqnonlin works much more reliable than fzero
	p  = lsqnonlin(@(p) spectral_density_bandpass_continuous_scale(fc,abs(p),pp) - Sc, p0);
%	p  = fzero(@(p) log(spectral_density_bandpass_continuous_scale(fc,abs(p),pp)) - log(Sc), p0);
	p  = abs(p);
end

