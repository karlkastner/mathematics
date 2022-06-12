% Mon 10 Jan 16:23:02 CET 2022
%
%% mode (maximum) of the lorentzian spectral density
%
function Sc = spectral_density_lorentzian_max(fc,p)
	IS = spectral_density_lorentzian_scale(fc,p);
	% because S has been normalized to 1 at fc
	Sc = 1./IS;
end
