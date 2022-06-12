% Mon 10 Jan 16:24:35 CET 2022
%
%% transform maximum of the lorentzian spectral density to its distribution parameters
%
function p = spectral_density_lorentzian_max2par(fc,Sc)
	p0 = 1;
	p = zeros(size(Sc));
	for idx=1:length(Sc)
		p(idx) = fzero(@(p) spectral_density_lorentzian_max(fc,p) - Sc(idx),p0);
	end
	p = abs(p);
end
