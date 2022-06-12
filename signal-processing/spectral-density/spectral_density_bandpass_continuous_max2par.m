% Fri  7 Jan 16:11:57 CET 2022
%
%% transform mode (maxima) of the bandpass spectral density into the paramter
%% of the underlying distribution 
%
% function [p] = spectral_density_bandpass_max2par(fc,Sc,p0)
function [p] = spectral_density_bandpass_continuous_max2par(fc,Sc,p0,q)
	if (nargin()<3)
		p0 = 1;
	end
	% log for better conditioning
	% the max is quite sensitive to small changes in p for large p
	if (0)
	lp  = fzero(@(p) spectral_density_bandpass_continuous_max(fc,exp(abs(p))) - Sc, log(p0));
	p = exp(lp);
	p = abs(p);
	else
	p  = fzero(@(p) spectral_density_bandpass_continuous_max(fc,abs(p),q) - Sc, p0);
	p  = abs(p);
	end
end

