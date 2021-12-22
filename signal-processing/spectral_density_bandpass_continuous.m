% Tue 23 Nov 17:31:17 CET 2021
% Karl Kastner, Berlin
%
% S : spectral density of the bandpass filter in continuos space
%     limit case of the discrete bandpass for dx -> 0
%
% f     : frequency (abszissa)
% f0    : central frequncy, location of maximum on abszissa
% order : number of times filter is applied iteratively, not necessarily integer
function S = spectral_density_bandpass_continous(fx,f0,order)
	S = (abs(fx)./(fx.^2 + f0.^2)).^(2*order);
	% normalize
	if (issym(fx))
		% note : this is only the correct scaling when p = 1
		S = 4/pi*f0*S;
	else
		% normalize
		I = spectral_density_area(fx,S);
		S = S./I;
	end
end

